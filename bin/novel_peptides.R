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
#'
#' Output:
#'   - {sample_name}_novel_peptides.tsv with columns: Sequence, PSM, Intensity,
#'     rna_detection_status, peptide_status, transcript_id, gene_id,
#'     n_transcripts, n_high_confidence, fasta_headers
#'     
#' Usage:
#'   Rscript novel_peptide.R --ms_search_software fragpipe --sample_name A549 --acquisition_type DIA --outdir results
#'   


# =============================================================================
# Load libraries
# =============================================================================

suppressPackageStartupMessages({
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

# =============================================================================
# Read inputs
# =============================================================================

# Fragpipe DIA
if (opt$acquisition_type == "DIA" & opt$ms_search_software == "fragpipe") {
  
  # peptides
  peptides_path = file.path("S5_PROTEOMICS/M2_FRAGPIPE", opt$sample_name, "peptide.tsv")
  dia_peptides = read_tsv(peptides_path, show_col_types = FALSE)
  
  peptides_df = dia_peptides %>%
    select(
      Sequence = `Peptide`,                   # Clean peptide sequence
      PSM = contains("Spectral Count"),       # PSM count
      Protein, 
      Additional_Proteins = `Mapped Proteins` # All protein matches (comma delimiter-separated)
    )
  
  # get intensities
  dia_path = file.path("S5_PROTEOMICS/M2_FRAGPIPE", opt$sample_name, "dia-quant-output", "report.tsv")
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
  
  peptides_path = file.path("S5_PROTEOMICS/M2_FRAGPIPE", opt$sample_name, "combined_peptide.tsv")
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
    n_transcripts = n(),
    n_high_confidence = sum(status == "high_confidence", na.rm = TRUE),
    fasta_headers = paste0(header, collapse = ",")
  ) %>%
  ungroup() %>%
  mutate(n_high_confidence = if_else(rna_detection_status == "RNA_not_detected", NA_real_, n_high_confidence)) %>%
  select(-all_proteins)

# =============================================================================
# Write summary and output
# =============================================================================

# Write output
outfile = file.path(opt$outdir, paste0(opt$sample_name, "_novel_peptides.tsv"))
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

cat("\nDone.\n")