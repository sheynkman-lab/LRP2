#!/usr/bin/env Rscript

#' Classifies peptides from mass spec database search results (Fragpipe)
#' 
#' Classification logic:
#'   - If a peptide maps to both LR and GENCODE proteins, GENCODE entries are
#'     discarded and only LR transcript IDs are retained. The peptide is
#'     classified as "annotated" (confirmed by both sources).
#'   - If a peptide maps only to LR proteins, it is classified as "novel"
#'     and all LR transcript IDs are retained.
#'   - If a peptide maps only to GENCODE proteins, it is classified as
#'     "annotated" with rna_detection_status = "RNA_not_detected". GENCODE
#'     transcript IDs are retained.
#'
#' Inputs:
#'   - Peptide search results from FragPipe (DDA or DIA) or MetaMorpheus
#'   - Optional: GTF with CDS coords for peptide mapping to genome
#'
#' Output:
#'   - {sample_name}.proteomics.novel_peptides.tsv with columns: Sequence, PSM, Intensity,
#'     rna_detection_status, peptide_status, transcript_id, gene_id,
#'     n_transcripts, n_high_confidence, fasta_headers
#'   - {sample_name}.proteomics.all_peptides.bed (optional)
#'     
#' Usage:
#'   Rscript novel_peptide.R --ms_search_software fragpipe --sample_name A549 --acquisition_type DIA --outdir results
#'   


# =============================================================================
# Load libraries
# =============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(tidyverse)
  library(magrittr)
  library(optparse)
})

options(scipen = 999)


# =============================================================================
# Command Line Arguments
# =============================================================================

option_list = list(
  make_option(c("--sample_name"), type = "character", default = NULL,
              help = "Sample name, used to locate search output [required]"),
  make_option(c("--ms_search_software"), type = "character", default = "fragpipe",
              help = "Search engine used: 'fragpipe' or 'metamorpheus' [default: fragpipe]"),
  make_option(c("--acquisition_type"), type = "character", default = NULL,
              help = "Fragpipe acquisition tpye: 'DIA' or 'DDA' [required for fragpipe]"),
  make_option(c("--fragpipe_results_dir"), type = "character", default = NULL,
              help = "Path to FragPipe results directory (e.g., results/S5_PROTEOMICS/M3_FRAGPIPE/sample_name)"),
  make_option(c("--custom_gtf"), type = "character", default = NULL,
              help = "Path to custom/LRP transcript GTF with CDS entries"),
  make_option(c("--custom_fasta"), type = "character", default = NULL,
              help = "Path to custom/LRP protein FASTA (headers: transcript_id|gene_id at minimum)"),
  make_option(c("--gencode_fasta"), type = "character", default = NULL,
              help = "Path to GENCODE protein FASTA (pc_translations.fa)"),
  make_option(c("--gencode_gtf"), type = "character", default = NULL,
              help = "Path to GENCODE annotation GTF (basic.annotation.gtf)"),
  make_option(c("--gencode_version"), type = "character", default = "46",
              help = "GENCODE version for GTF download [default: 46]"),
  make_option(c("--species"), type = "character", default = "human",
              help = "Species: 'human' or 'mouse' [default: human]"),
  make_option(c("--outdir"), type = "character", default = ".",
              help = "Output directory [default: .]")
)

opt = parse_args(OptionParser(option_list = option_list))

# =============================================================================
# Validate inputs
# =============================================================================

if (is.null(opt$sample_name)) stop("--sample_name is required")
if (!opt$ms_search_software %in% c("fragpipe", "metamorpheus")) {
  stop("--ms_search_software must be 'fragpipe' or 'metamorpheus'")
}

# Validate acquisition_type is provided when using fragpipe
if (opt$ms_search_software == "fragpipe" & is.null(opt$acquisition_type)) {
  stop("--acquisition_type is required when using fragpipe (must be 'DDA' or 'DIA')")
}
if (!is.null(opt$acquisition_type) && !opt$acquisition_type %in% c("DDA", "DIA")) {
  stop("--acquisition_type must be 'DDA' or 'DIA'")
}

if (!opt$species %in% c("human", "mouse")) {
  stop("--species must be 'human' or 'mouse' to pull proper GENCODE gtf")
}

# optional custom fasta and gtf requirements if peptide mapping
if (!is.null(opt$custom_gtf) && is.null(opt$custom_fasta)) {
  stop("--custom_fasta is required when --custom_gtf is provided")
}
if (!is.null(opt$custom_fasta) && is.null(opt$custom_gtf)) {
  stop("--custom_gtf is required when --custom_fasta is provided")
}

# gencode fasta and gtf should be provided together
if (!is.null(opt$gencode_fasta) && is.null(opt$gencode_gtf)) {
  stop("--gencode_gtf is required when --gencode_fasta is provided")
}
if (!is.null(opt$gencode_gtf) && is.null(opt$gencode_fasta)) {
  stop("--gencode_fasta is required when --gencode_gtf is provided")
}

# =============================================================================
# Function to map peptides to genome
# =============================================================================

#' Map peptides to genomic coordinates (BED12)
#'
#' Takes peptide sequences with transcript IDs, locates each peptide within
#' the corresponding protein sequence, then maps to genomic coordinates
#' using CDS blocks from a GTF. Output is BED12 format.
#'
#' Peptides spanning exon-exon junctions are represented as multi-block
#' BED12 entries.
#'
#' @param peptides Tibble with at minimum Sequence and transcript_id columns.
#'   Only distinct Sequence + transcript_id pairs are used.
#' @param aa_fasta_path Path to amino acid protein FASTA or AAStringSet object.
#' @param gtf Path to GTF with CDS entries. Must have columns: seqnames, start,
#'   end, strand, type, and transcript_id.
#' @return Tibble in BED12 format with extra columns: Sequence, transcript_id,
#'   aa_start, aa_end.
#'
#' @examples
#' bed = map_peptides_to_genome(peptides, "gencode.v46.pc_translations.fa", gtf)
#' write_tsv(bed %>% select(1:12), "peptides.bed", col_names = FALSE)
map_peptides_to_genome = function(peptides, aa_fasta_path, gtf_path) {
  
  # --- 0. Input validation ---
  
  # Check peptide table structure
  peptide_tx = peptides 
  if (ncol(peptide_tx) < 2) stop("Peptide table must have at least 2 columns (Sequence, transcript_id)")
  if (!is.character(peptide_tx[[1]])) stop("First column of peptides table must be character (peptide sequences)")
  if (!is.character(peptide_tx[[2]])) stop("Second column of peptides table must be character (transcript IDs)")
  names(peptide_tx)[1:2] = c("Sequence", "transcript_id")
  has_psm = ncol(peptide_tx) >= 3 && is.numeric(peptide_tx[[3]])
  if (has_psm) names(peptide_tx)[3] = "PSM"
  
  # Check GTF has CDS entries
  gtf = import(gtf_path) %>% as.data.frame()
  if (!"type" %in% names(gtf)) stop("GTF must have a 'type' column")
  if (!any(gtf$type == "CDS")) stop("GTF has no CDS entries - cannot map peptides to genomic coordinates")
  
  # Check GTF has required columns
  required_gtf_cols = c("seqnames", "start", "end", "strand", "transcript_id")
  missing_gtf = setdiff(required_gtf_cols, names(gtf))
  if (length(missing_gtf) > 0) stop("GTF missing required columns: ", paste(missing_gtf, collapse = ", "))
  
  # --- 1. Read protein sequences ---
  aa_fasta = readAAStringSet(aa_fasta_path)

  # Detect FASTA format from first header
  # GENCODE format: ENSP|ENST|ENSG|... (protein_id|transcript_id|gene_id)
  # LR format: GENE::TRANSCRIPT|GENE_ID|... (transcript_id|gene_id)
  first_header = names(aa_fasta)[1]
  is_gencode_format = str_detect(first_header, "^ENSP")

  if (is_gencode_format) {
    # GENCODE: extract second field (transcript_id)
    proteins = tibble(
      header = names(aa_fasta),
      aa_seq = as.character(aa_fasta)
    ) %>%
      separate_wider_delim(header, delim = "|", names = c("protein_id", "transcript_id", "gene_id"),
                           too_many = "drop", cols_remove = FALSE) %>%
      select(-protein_id, -gene_id)  # Keep only transcript_id
  } else {
    # LR: extract first field (transcript_id)
    proteins = tibble(
      header = names(aa_fasta),
      aa_seq = as.character(aa_fasta)
    ) %>%
      separate_wider_delim(header, delim = "|", names = c("transcript_id", "gene_id"),
                           too_many = "drop", cols_remove = FALSE) %>%
      select(-gene_id)  # Keep transcript_id, drop gene_id
  }
  
  # Check transcript ID overlap across inputs
  tx_in_peptides = unique(peptide_tx$transcript_id)
  tx_in_fasta    = unique(proteins$transcript_id)
  tx_in_gtf      = unique(gtf$transcript_id[gtf$type == "CDS"])
  
  overlap_fasta    = sum(tx_in_peptides %in% tx_in_fasta)
  overlap_gtf      = sum(tx_in_peptides %in% tx_in_gtf)
  fasta_not_in_gtf = sum(!tx_in_fasta %in% tx_in_gtf)
  gtf_not_in_fasta = sum(!tx_in_gtf %in% tx_in_fasta)
  
  cat("Transcript ID overlap:\n")
  cat("  Peptide Transcripts in FASTA:", overlap_fasta, "/", length(tx_in_peptides), "\n")
  cat("  Peptide Transcripts in CDS GTF:", overlap_gtf, "/", length(tx_in_peptides), "\n")

  if (fasta_not_in_gtf > 0) cat("  Warning:", fasta_not_in_gtf, "FASTA transcripts missing from GTF CDS\n")
  if (gtf_not_in_fasta > 0) cat("  Warning:", gtf_not_in_fasta, "GTF CDS transcripts missing from FASTA\n")
  if (overlap_fasta == 0) stop("No transcript IDs in peptide table match the FASTA")
  if (overlap_gtf == 0) stop("No transcript IDs in peptide table match CDS entries in the GTF")
  
  # --- 2. Locate peptide within protein sequence ---
  peptide_positions = peptide_tx %>%
    left_join(proteins, by = "transcript_id") %>%
    mutate(
      aa_start = str_locate(aa_seq, fixed(Sequence))[, 1],
      aa_end   = str_locate(aa_seq, fixed(Sequence))[, 2]
    ) %>%
    filter(!is.na(aa_start)) %>%
    mutate(
      cds_nt_start = (aa_start - 1) * 3 + 1,
      cds_nt_end   = aa_end * 3
    ) %>%
    select(-aa_seq, -header)

  # --- 3. Build CDS cumulative lengths from GTF ---
  cds_blocks = gtf %>%
    filter(type == "CDS") %>%
    select(seqnames, start, end, strand, transcript_id) %>%
    arrange(transcript_id, if_else(strand == "+", start, -start)) %>%
    mutate(cds_length = end - start + 1) %>%
    group_by(transcript_id) %>%
    mutate(
      cumulative_length = cumsum(cds_length),
      prior_length = cumulative_length - cds_length
    ) %>%
    ungroup()
  
  # --- 4. Find all CDS blocks that overlap the peptide ---
  peptide_blocks = peptide_positions %>%
    left_join(cds_blocks, by = "transcript_id", relationship = "many-to-many") %>%
    filter(cds_nt_end > prior_length & cds_nt_start <= cumulative_length) %>%
    mutate(
      # clip block to peptide boundaries
      block_start_nt = pmax(cds_nt_start, prior_length + 1),
      block_end_nt   = pmin(cds_nt_end, cumulative_length),
      
      # convert to genomic coordinates
      genomic_block_start = if_else(strand == "+",
                                    start + (block_start_nt - prior_length) - 1,
                                    end - (block_end_nt - prior_length) + 1),
      genomic_block_end = if_else(strand == "+",
                                  start + (block_end_nt - prior_length) - 1,
                                  end - (block_start_nt - prior_length) + 1)
    )
  
  # --- 5. Collapse to BED12 format ---
  # BED12 requires blocks sorted by genomic position ascending
  bed12 = peptide_blocks %>%
    arrange(Sequence, transcript_id, genomic_block_start) %>%
    group_by(Sequence, transcript_id, seqnames, strand, aa_start, aa_end) %>%
    summarize(
      chrom_start  = min(genomic_block_start) - 1,  # BED is 0-based half-open
      chrom_end    = max(genomic_block_end),
      block_count  = n(),
      block_sizes  = paste0(genomic_block_end - genomic_block_start + 1, collapse = ","),
      block_starts = paste0(genomic_block_start - min(genomic_block_start), collapse = ","),
      .groups = "drop"
    ) %>%
    transmute(
      chrom       = seqnames,
      chrom_start,
      chrom_end,
      name        = Sequence,
      score       = 0,
      strand,
      thick_start = chrom_start,
      thick_end   = chrom_end,
      item_rgb    = "0,0,0",
      block_count,
      block_sizes,
      block_starts,
      Sequence,
      transcript_id,
      aa_start,
      aa_end
    )
  
  # Join PSM back if provided
  if (has_psm) {
    psm_vals = peptide_tx %>%
      distinct(Sequence, transcript_id, PSM)
    bed12 = bed12 %>%
      left_join(psm_vals, by = c("Sequence", "transcript_id")) %>%
      mutate(name = paste0(Sequence, "|PSM=", PSM))
  }
  
  cat("\nMapped", nrow(bed12), "peptides to genomic coordinates\n")
  n_junction = sum(bed12$block_count > 1)
  if (n_junction > 0) cat("  ", n_junction, "span exon-exon junctions\n")
  
  unmapped = peptide_tx %>%
    anti_join(bed12, by = c("Sequence", "transcript_id"))
  
  if (nrow(unmapped) > 0) {
    cat("  Warning:", nrow(unmapped), "peptides could not be mapped to genomic coordinates\n")
  }
  
  return(bed12)
}

# =============================================================================
# Read inputs
# =============================================================================

# Fragpipe DIA
if (opt$acquisition_type == "DIA" & opt$ms_search_software == "fragpipe") {

  # Determine base path for FragPipe results
  if (!is.null(opt$fragpipe_results_dir)) {
    fragpipe_base = opt$fragpipe_results_dir
  } else {
    fragpipe_base = file.path("S5_PROTEOMICS/M3_FRAGPIPE", opt$sample_name)
  }

  # peptides
  peptides_path = file.path(fragpipe_base, "peptide.tsv")
  dia_peptides = read_tsv(peptides_path, show_col_types = FALSE)

  peptides_df = dia_peptides %>%
    select(
      Sequence = `Peptide`,                   # Clean peptide sequence
      PSM = contains("Spectral Count"),       # PSM count
      Protein,
      Additional_Proteins = `Mapped Proteins` # All protein matches (comma delimiter-separated)
    )

  # get intensities
  dia_path = file.path(fragpipe_base, "dia-quant-output", "report.tsv")
  dia_report = read_tsv(dia_path, show_col_types = FALSE)
  
  dia_intensity = dia_report %>%
    select(Sequence = Stripped.Sequence, Precursor.Quantity) %>%
    group_by(Sequence) %>%
    summarize(Intensity = sum(Precursor.Quantity, na.rm = TRUE), .groups = "drop")
  
  peptides_df %<>%
    left_join(dia_intensity, by = c("Sequence")) %>%
    mutate(Intensity = replace_na(Intensity, 0))
  
}

# Fragpipe DDA
if (opt$acquisition_type == "DDA" & opt$ms_search_software == "fragpipe") {

  # Determine base path for FragPipe results
  if (!is.null(opt$fragpipe_results_dir)) {
    fragpipe_base = opt$fragpipe_results_dir
  } else {
    fragpipe_base = file.path("S5_PROTEOMICS/M3_FRAGPIPE", opt$sample_name)
  }

  peptides_path = file.path(fragpipe_base, "combined_peptide.tsv")
  dda_peptides = read_tsv(peptides_path, show_col_types = FALSE)
  
  peptides_df = dda_peptides %>%
    select(
      Sequence = `Peptide Sequence`,          # Clean peptide sequence
      PSM = contains("Spectral Count"),       # PSM count
      Intensity = contains("Intensity"),      # Ion intensity
      Protein, 
      Additional_Proteins = `Mapped Proteins` # All protein matches (comma delimiter-separated)
    )
  
}

# =============================================================================
# Identify novel peptides
# =============================================================================

peptides_longer = peptides_df %>% 
  mutate(all_proteins = if_else(
    Additional_Proteins == "" | is.na(Additional_Proteins),
    Protein,
    paste(Protein, Additional_Proteins, sep = ","))) %>%
  select(-Protein, -Additional_Proteins) %>%
  mutate(header = all_proteins) %>%
  separate_rows(header, sep = ",") %>%
  mutate(header = str_trim(header))

# logic to call novel peptides
status = peptides_longer %>%
  separate_wider_delim(header, 
                       delim = "|",
                       names = c("transcript_id", "gene_id", "pclass", "status", "reference_type"), 
                       cols_remove = FALSE) %>%
  group_by(Sequence, all_proteins) %>%
  mutate(
    maps_any_lrs = any(reference_type != "gencode", na.rm = TRUE),
    maps_any_ref = any(reference_type == "gencode", na.rm = TRUE)
  ) %>%
  
  # discard gencode rows when LR is also present
  filter(!(maps_any_lrs & reference_type == "gencode")) %>%
  mutate(
    peptide_status = case_when(
      maps_any_lrs & maps_any_ref ~ "annotated",
      maps_any_lrs & !maps_any_ref ~ "novel",
      !maps_any_lrs & maps_any_ref ~ "annotated"
    ),
    rna_detection_status = case_when(
      maps_any_lrs ~ "RNA_detected",
      !maps_any_lrs & maps_any_ref ~ "RNA_not_detected"
    )
  ) %>%
  ungroup()

out = status %>%
  group_by(Sequence, all_proteins, PSM, Intensity, rna_detection_status, peptide_status) %>%
  arrange(pclass) %>%
  summarize(
    transcript_id = paste0(transcript_id, collapse = ","),
    gene_id = paste0(unique(gene_id), collapse = ","),
    best_status = case_when(
      any(status == "high_confidence") ~ "high_confidence",
      any(status == "sqanti_atypical") ~ "sqanti_atypical",
      any(status == "NMD") ~ "NMD"
    ),
    n_transcripts = n(),
    n_high_confidence = sum(status == "high_confidence", na.rm = TRUE),
    fasta_headers = paste0(header, collapse = ",")
  ) %>%
  ungroup() %>%
  mutate(n_high_confidence = if_else(rna_detection_status == "RNA_not_detected", NA_real_, n_high_confidence),
         best_status = if_else(rna_detection_status == "RNA_not_detected", NA_character_, best_status)) %>%
  select(-all_proteins)

# =============================================================================
# Write novel peptide summary and output 
# =============================================================================

# Write output
outfile = file.path(opt$outdir, paste0(opt$sample_name, ".proteomics.novel_peptides.tsv"))
write_tsv(out, outfile)

# Summary
cat("\n=== Novel Peptide Summary ===\n")
cat("Total peptides:", nrow(out), "\n")

n_rna_detected = sum(out$rna_detection_status == "RNA_detected", na.rm = TRUE)
n_rna_not_detected = sum(out$rna_detection_status == "RNA_not_detected", na.rm = TRUE)
n_novel = sum(out$peptide_status == "novel", na.rm = TRUE)
n_annotated_rna = sum(out$peptide_status == "annotated" & out$rna_detection_status == "RNA_detected", na.rm = TRUE)
n_annotated_gencode = sum(out$peptide_status == "annotated" & out$rna_detection_status == "RNA_not_detected", na.rm = TRUE)

cat("\nRNA detection:\n")
cat("  Detected at RNA level:", n_rna_detected, "\n")
cat("  GENCODE only (no RNA):", n_rna_not_detected, "\n")

cat("\nPeptide classification:\n")
cat("  Novel (LR only):", n_novel, "\n")
cat("  Annotated (LR + GENCODE):", n_annotated_rna, "\n")
cat("  Annotated (GENCODE only):", n_annotated_gencode, "\n")
cat("  Peptides with best status NMD:", sum(out$best_status == "NMD", na.rm = TRUE), "\n")

# =============================================================================
# Map peptides to genomic coordinates
# =============================================================================

# GENCODE and custom/LRP are independent — either or both can run
# Each requires its own FASTA + GTF pair

gencode_bed = NULL
custom_bed = NULL

# --- GENCODE peptides (runs whenever there are GENCODE-only peptides) ---

gencode_peptides = out %>%
  filter(rna_detection_status == "RNA_not_detected")

if (nrow(gencode_peptides) > 0) {
  cat("\n=== GENCODE peptides ===\n")
  
  gencode_peptides %<>%
    select(Sequence, transcript_id, gene_id, PSM) %>%
    separate_rows(transcript_id, sep = ",") %>%
    separate_rows(gene_id, sep = ",") %>%
    filter(transcript_id != "NA")
  
  # Remove peptides mapping to multiple genes
  multi_gene = gencode_peptides %>%
    distinct(Sequence, gene_id) %>%
    count(Sequence) %>%
    filter(n > 1)
  cat("Removing", nrow(multi_gene), "peptides mapping to multiple genes\n")
  
  unique_gencode_peptides = gencode_peptides %>%
    filter(!Sequence %in% multi_gene$Sequence) %>%
    distinct(Sequence, transcript_id, PSM) %>%
    group_by(Sequence, PSM) %>%
    slice(1) %>%
    ungroup()
    
  # Use provided GENCODE files, fall back to download
  if (!is.null(opt$gencode_fasta) && !is.null(opt$gencode_gtf)) {
    gencode_fa_local_path = opt$gencode_fasta
    gencode_gtf_local_path = opt$gencode_gtf
  } else {
    if (opt$species == "human") {
      gencode_base = paste0("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_", opt$gencode_version)
      gencode_ver = paste0("v", opt$gencode_version)
    } else {
      gencode_base = paste0("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M", opt$gencode_version)
      gencode_ver = paste0("vM", opt$gencode_version)
    }
    
    gencode_fasta = file.path(gencode_base, paste0("gencode.", gencode_ver, ".pc_translations.fa.gz"))
    gencode_gtf = file.path(gencode_base, paste0("gencode.", gencode_ver, ".basic.annotation.gtf.gz"))
    
    gencode_fa_local_path = file.path(opt$outdir, basename(gencode_fasta))
    gencode_gtf_local_path = file.path(opt$outdir, basename(gencode_gtf))
    
    if (!file.exists(gencode_fa_local_path)) {
      download.file(gencode_fasta, gencode_fa_local_path, mode = "wb")
    }
    if (!file.exists(gencode_gtf_local_path)) {
      download.file(gencode_gtf, gencode_gtf_local_path, mode = "wb")
    }
  }
  
  # Map GENCODE peptides
  gencode_bed = map_peptides_to_genome(unique_gencode_peptides, gencode_fa_local_path, gencode_gtf_local_path)
  
}
  
# --- Custom/LRP peptides (only when custom/lrp fasta + gtf provided) ---
if (!is.null(opt$custom_fasta) && !is.null(opt$custom_gtf)) {
  cat("\n=== Custom/LRP Peptides ===\n")
  
  lr_peptides = out %>%
    filter(rna_detection_status == "RNA_detected") %>%
    select(Sequence, transcript_id, gene_id, PSM) %>%
    separate_rows(transcript_id, sep = ",") %>%
    separate_rows(gene_id, sep = ",") %>%
    filter(transcript_id != "NA")
  
  # Remove peptides mapping to multiple genes
  multi_gene_lr = lr_peptides %>%
    distinct(Sequence, gene_id) %>%
    count(Sequence) %>%
    filter(n > 1)
  cat("\nRemoving", nrow(multi_gene_lr), "LR peptides mapping to multiple genes\n")
  
  unique_lr_peptides = lr_peptides %>%
    filter(!Sequence %in% multi_gene_lr$Sequence) %>%
    distinct(Sequence, transcript_id, PSM) %>%
    group_by(Sequence, PSM) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(transcript_id = as.character(transcript_id))
  
  # Map custom/LRP peptides
  custom_bed = map_peptides_to_genome(unique_lr_peptides, opt$custom_fasta, opt$custom_gtf)
}

all_bed = bind_rows(gencode_bed, custom_bed)

if (nrow(all_bed) > 0) {
  bed_outfile = file.path(opt$outdir, paste0(opt$sample_name, ".proteomics.all_peptides.bed"))
  write_tsv(all_bed %>% select(1:12), bed_outfile, col_names = FALSE)
  cat("\nWrote", nrow(all_bed), "BED12 entries to", bed_outfile, "\n")
}


cat("\nDone.\n")