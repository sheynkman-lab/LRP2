#!/usr/bin/env Rscript

#' SQANTI3 Transcript Filtering Script
#' 
#' Generate hash ids for all transcripts before filtering
#' 
#' R script to filter SQANTI3 transcripts based on quality criteria:
#' - Protein-coding genes only
#' - Remove transcripts with internal priming 
#' - Remove template switching artifacts
#' - Filter by SQANTI structural categories
#'
#' Outputs both filtered files and dropout files with reasons.

# =============================================================================
# Load required libraries
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
  library(rtracklayer)
  library(Biostrings) 
  library(data.table)
  library(magrittr)
})

# =============================================================================
# Get environment variables for filtering parameters
# =============================================================================

dir                     <- Sys.getenv("OUTPUT_DIR")
basename_prefix         <- Sys.getenv("OUTPUT_BASE_NAME")
basename                <- paste0(basename_prefix, "_sqanti_transcript")
gencode_gtf             <- Sys.getenv("GENCODE_GTF_FILE")
metadata_path           <- Sys.getenv("SAMPLE_METADATA")
filter_protein_coding   <- as.logical(Sys.getenv("PROTEIN_CODING_FILTER", "TRUE"))
filter_internal_priming <- as.logical(Sys.getenv("INTERNAL_PRIMING_FILTER", "TRUE"))
filter_RTS              <- as.logical(Sys.getenv("TEMPLATE_SWITCHING_FILTER", "TRUE"))
percent_polyA_threshold <- as.numeric(Sys.getenv("PERCENT_POLYA_THRESHOLD", "60"))
structural_level        <- Sys.getenv("STRUCTURE_FILTER", "strict")

# =============================================================================
# Check required files
# =============================================================================

classification_file = paste0(dir, "/", basename, "_classification.txt")
sqanti_gtf          = paste0(dir, "/", basename, "_corrected.gtf")
sqanti_fasta        = paste0(dir, "/", basename, "_corrected.fasta")

stopifnot("File not found" = file.exists(classification_file))
stopifnot("File not found" = file.exists(sqanti_gtf))
stopifnot("File not found" = file.exists(sqanti_fasta))

# Create output directory
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
dropout_dir = dir

# =============================================================================
# Helper functions
# =============================================================================

#' Convert gtf to psl format- adapted from FLAIR and IsoViz
#' @param gtf_input_path Input GTF file path
#' @param psl_output_file Output PSL file
convert_gtf_to_psl = function(gtf_input_path, psl_output_file){
  
  # Read the GTF file 
  gr     = import(gtf_input_path, format = "gtf")
  gtf_df = as.data.frame(gr)
  
  # Set options to avoid scientific notation
  options(scipen = 999)
  
  # keep only the following columns 
  gtf_df %<>% 
    dplyr::select(seqnames, type, start, end, strand, transcript_id, gene_id)
  
  colnames(gtf_df) = c("chrom", "ty", "start", "end", "strand", "transcript_id", "gene_id")
  gtf_df$start     = gtf_df$start-1
  
  # Filter for exons only
  exons_df = gtf_df %>% 
    filter(ty == "exon")
  
  cols_to_convert           = c("chrom", "ty", "strand", "transcript_id", "gene_id")
  exons_df[cols_to_convert] = lapply(exons_df[cols_to_convert], as.character)
  
  # Function to generate PSL lines from grouped exons
  generate_psl_lines = function(exons_grouped) {
    
    psl_lines = lapply(exons_grouped, function(test_ex) {
      
      # Get block starts and sizes
      blockstarts = test_ex$start
      blocksizes  = test_ex$end - test_ex$start
      blockcount  = length(blockstarts)
      
      # Reverse if needed (negative strand)
      if (blockcount > 1 && blockstarts[1] > blockstarts[2]) {
        blocksizes  = rev(blocksizes)
        blockstarts = rev(blockstarts)
      }
      
      # Calculate transcript coordinates
      tstart = blockstarts[1]
      tend   = blockstarts[blockcount] + blocksizes[blockcount]
      qsize  = sum(blocksizes)
      qname  = paste0(test_ex$transcript_id[1], "_", test_ex$gene_id[1])
      
      # Calculate query starts (cumulative positions)
      qstarts = c(0, cumsum(blocksizes)[-blockcount])
      
      # Format as comma-separated strings
      qstarts_str     = paste0(paste(qstarts, collapse = ","), ",")
      blocksizes_str  = paste0(paste(blocksizes, collapse = ","), ",")
      blockstarts_str = paste0(paste(blockstarts, collapse = ","), ",")
      
      # Construct PSL line
      psl_line = c(
        0, 0, 0, 0, 0, 0, 0, 0, 
        test_ex$strand[1], qname, qsize, 0, qsize,
        test_ex$chrom[1], 0, tstart, tend, blockcount, 
        blocksizes_str, qstarts_str, blockstarts_str
      )
      
      # Return as tab-separated string
      paste(psl_line, collapse = "\t")
    })
    
    # Return all lines
    unlist(psl_lines)
  }
  
  # Usage:
  exons_grouped = exons_df %>%
    group_by(transcript_id) %>% 
    group_split()
  
  psl_lines = generate_psl_lines(exons_grouped)
  
  writeLines(psl_lines, psl_output_file)
}

#' Save filtered sequences to FASTA
#' @param fasta_file Input FASTA file
#' @param keep_ids IDs to keep
#' @param output_file Output FASTA file
save_filtered_fasta <- function(fasta_file, keep_ids, output_file) {

  sequences          = readDNAStringSet(fasta_file)
  keep_mask          = names(sequences) %in% keep_ids
  filtered_sequences = sequences[keep_mask]
  
  # Write filtered FASTA
  writeXStringSet(filtered_sequences, output_file)
  message(paste0("Saved ", length(filtered_sequences), " sequences to ", basename(output_file)))
}

#' Save filtered GTF
#' @param gtf_file Input GTF file
#' @param keep_ids Transcript IDs to keep
#' @param output_file Output GTF file
save_filtered_gtf <- function(gtf_file, keep_ids, output_file) {

  gtf          = import(gtf_file)
  filtered_gtf = gtf[gtf$transcript_id %in% keep_ids]
  
  # Write filtered GTF
  export(filtered_gtf, output_file)
  message(paste0("Saved ", length(unique(filtered_gtf$transcript_id)), " transcripts to ", basename(output_file)))
}

# =============================================================================
# Generating hashids
# =============================================================================

message("\n Generating hash ids for all transcripts...")

# Sample gtf to psl
sqanti_psl = file.path(dir, paste0(basename, "_corrected.psl")) # output psl
convert_gtf_to_psl(gtf_input_path  = sqanti_gtf, 
                   psl_output_file = sqanti_psl)

# Run python script for hash ids- generates mapping file with transcript_id and hash_id
mapping_file = file.path(dir, paste0(basename, "_hashids_mapping.txt"))
# Set path to the hashing function Python script from environment for proper usage in module container
python_script = Sys.getenv("HASHLIB_SCRIPT", "hashlib_id_generator.py")
system2("python", args = c(python_script, sqanti_psl, mapping_file))

# =============================================================================
# Combine sqanti associated gene_ids with gencode ids and hash ids, include counts
# =============================================================================

# Gencode gtf used to select protein coding genes and add column of readable gene name
# Use the environment variable passed from Nextflow, not hardcoded path
gencode_gtf = Sys.getenv("GENCODE_GTF_FILE")
if (gencode_gtf == "") {
  stop("GENCODE_GTF_FILE environment variable not set")
}
gencode    = import(gencode_gtf, format = "gtf")
gencode_df = as.data.frame(gencode)

# Filter for transcript entries
gencode_df %<>% filter(type == "transcript")

# Rename columns based on what's available, handling gene_type/gene_biotype variations
gencode_df %<>%
  rename(
    associated_gene = gene_id,
    associated_transcript = transcript_id
  )

# Handle gene_type vs gene_biotype naming
if ("gene_biotype" %in% colnames(gencode_df) && !("gene_type" %in% colnames(gencode_df))) {
  gencode_df %<>% rename(gene_type = gene_biotype)
}

# Select only the columns we need, using only those that exist
select_cols <- c("associated_gene", "associated_transcript")
if ("gene_type" %in% colnames(gencode_df)) select_cols <- c(select_cols, "gene_type")
if ("gene_name" %in% colnames(gencode_df)) select_cols <- c(select_cols, "gene_name")
if ("transcript_name" %in% colnames(gencode_df)) select_cols <- c(select_cols, "transcript_name")
if ("protein_id" %in% colnames(gencode_df)) select_cols <- c(select_cols, "protein_id")

gencode_df %<>% select(all_of(select_cols))

# Select columns for gencode_gene, handling missing columns
gencode_gene_cols <- c("associated_gene")
if ("gene_type" %in% colnames(gencode_df)) gencode_gene_cols <- c(gencode_gene_cols, "gene_type")
if ("gene_name" %in% colnames(gencode_df)) gencode_gene_cols <- c(gencode_gene_cols, "gene_name")

gencode_gene = gencode_df %>%
  select(all_of(gencode_gene_cols)) %>%
  distinct()

# Metatable to add specific sample ids
metadata      = fread(metadata_path)
sample_lookup = setNames(metadata$name, metadata$bam_id)

# SQANTI classification file- add sample names, calculate cpm on all transcripts before filtering
classification_file = file.path(dir, paste0(basename, "_classification.txt"))
sqanti_df = read_tsv(classification_file, show_col_types = FALSE)

colnames(sqanti_df) = sub("^FL\\.", "", colnames(sqanti_df))
colnames(sqanti_df) = ifelse(colnames(sqanti_df) %in% names(sample_lookup),
                             paste0(sample_lookup[colnames(sqanti_df)], "_counts"),
                             colnames(sqanti_df))

counts_cols = grep("_counts$", colnames(sqanti_df), value = TRUE)
sqanti_df %<>%
  mutate(across(all_of(counts_cols), 
                ~ (.x / sum(.x)) * 1e6,
                .names = "{gsub('_counts', '_cpm', .col)}"))

# Keep only isoforms with ENS gene IDs (if available, otherwise we keep all)
ens_genes <- sum(str_starts(sqanti_df$associated_gene, "ENS"), na.rm = TRUE)
if (ens_genes > 0) {
  sqanti_df_full = sqanti_df %>%
    filter(!is.na(associated_gene), str_starts(associated_gene, "ENS"))
} else {
  sqanti_df_full = sqanti_df %>%
    filter(!is.na(associated_gene))
}

# Build join columns dynamically based on what exists in gencode_gene
join_cols_1 <- "associated_gene"

# Left join with gencode_gene (has at least associated_gene)
sqanti_df_full %<>%
  left_join(gencode_gene, by = join_cols_1)

# Build join columns for second join with gencode_df
join_cols_2 <- c("associated_gene", "associated_transcript")
if ("gene_type" %in% colnames(gencode_df)) {
  join_cols_2 <- c(join_cols_2, "gene_type")
}
if ("gene_name" %in% colnames(gencode_df)) {
  join_cols_2 <- c(join_cols_2, "gene_name")
}

# Second join only if we have columns beyond the basics
if (length(join_cols_2) > 2) {
  sqanti_df_full %<>%
    left_join(gencode_df, by = join_cols_2)
}

# combine with hashids, and selct cpm columns
hashids = read_tsv(mapping_file)

# Build select columns dynamically based on what exists in sqanti_df_full
select_cols_final <- c("isoform_id" = "isoform", "ensg_gene_id" = "associated_gene", "enst_transcript_id" = "associated_transcript")

if ("gene_name" %in% colnames(sqanti_df_full)) select_cols_final <- c(select_cols_final, "gene_name" = "gene_name")
if ("transcript_name" %in% colnames(sqanti_df_full)) select_cols_final <- c(select_cols_final, "transcript_name" = "transcript_name")
if ("gene_type" %in% colnames(sqanti_df_full)) select_cols_final <- c(select_cols_final, "gene_type" = "gene_type")
if ("protein_id" %in% colnames(sqanti_df_full)) select_cols_final <- c(select_cols_final, "protein_id" = "protein_id")

sqanti_ids = sqanti_df_full %>%
  select(all_of(select_cols_final), ends_with("_counts"), ends_with("_cpm"))

all_ids = hashids %>% select(isoform_id = transcript_id, hash_id) %>%
  inner_join(sqanti_ids, by = "isoform_id") 

write_tsv(all_ids, file.path(dir, paste0(basename, "_all_hashids_with_cpm.txt")))

# =============================================================================
# Apply filters sequentially
# =============================================================================

message("\n Starting filtering of SQANTI output...")

# Initialize dropout tracking
original_ids    = sqanti_df_full$isoform
dropout_tracker = list()

# Keep only transcripts of protein coding genes
sqanti_df = sqanti_df_full
if (filter_protein_coding) {
  if ("gene_type" %in% colnames(sqanti_df)) {
    ids = sqanti_df$isoform
    sqanti_df %<>%
      filter(gene_type == "protein_coding")

    kept_ids    = sqanti_df$isoform
    dropped_ids = setdiff(ids, kept_ids)
    dropout_tracker[["not_protein_coding"]] = dropped_ids

    message("Protein coding filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts")
  }
}

# Filter based on downstream polyA sequence to remove internal priming
if (filter_internal_priming) {
  ids = sqanti_df$isoform

  # keep all isoforms below threshold
  keep1 = sqanti_df %>%
    filter(perc_A_downstream_TTS <= percent_polyA_threshold)

  # Get downstream sequence for FSM above threshold
  protected_seqs = sqanti_df %>%
    filter(perc_A_downstream_TTS > percent_polyA_threshold, structural_category == "full-splice_match") %>%
    distinct(associated_gene, seq_A_downstream_TTS)

  # Keep high polyA isoforms that share sequence with FSM transcript
  keep2 = sqanti_df %>%
    inner_join(protected_seqs, by = c("associated_gene", "seq_A_downstream_TTS"))

  sqanti_df = bind_rows(keep1, keep2) %>%
    distinct()

  kept_ids    = sqanti_df$isoform
  dropped_ids = setdiff(ids, kept_ids)
  dropout_tracker[["internal_priming"]] = dropped_ids

  message("Internal priming filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts")
}

# Filter out template switching, keep RTS_stage = FALSE (no RTS)
if (filter_RTS) {
  ids = sqanti_df$isoform

  # Find column positions
  ratio_tss_pos = which(colnames(sqanti_df) == "ratio_TSS")
  gene_type_pos = which(colnames(sqanti_df) == "gene_type")

  # Get columns between them (if any exist) and sum
  if (length(ratio_tss_pos) > 0 && length(gene_type_pos) > 0 && gene_type_pos - ratio_tss_pos > 1) {
    count_cols           = (ratio_tss_pos + 1):(gene_type_pos - 1)
    sqanti_df$avg_counts = rowMeans(sqanti_df[, count_cols], na.rm = TRUE)
  } else {
    sqanti_df$avg_counts = NA  # No columns in between
  }

  # filtering considerations
  sqanti_df %<>%
    filter(
      structural_category == "full-splice_match" |
      RTS_stage == FALSE |
      (RTS_stage == TRUE & all_canonical == TRUE & avg_counts > 3)
      )

  kept_ids    = sqanti_df$isoform
  dropped_ids = setdiff(ids, kept_ids)
  dropout_tracker[["Template_switching_artifact"]] = dropped_ids

  message("Template switching filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts")
}

# Apply structural category filter
STRUCTURAL_CATEGORIES <- list(
  strict = c("novel_not_in_catalog", "novel_in_catalog",
             "incomplete-splice_match", "full-splice_match"),
  all    = c("antisense", "novel_not_in_catalog", "novel_in_catalog",
             "incomplete-splice_match", "full-splice_match", "genic",
             "intergenic", "fusion", "genic_intron")
)

allowed_categories = STRUCTURAL_CATEGORIES[[structural_level]]
ids                = sqanti_df$isoform
sqanti_df %<>%
  filter(structural_category %in% allowed_categories)

kept_ids    = sqanti_df$isoform
dropped_ids = setdiff(ids, kept_ids)
dropout_tracker[["structural_category_filtered"]] = dropped_ids
message("Structural category filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts")

# Calculate final results
all_dropout_ids = setdiff(original_ids, kept_ids)

message(paste0("Final: ", length(kept_ids), " kept, ", length(all_dropout_ids), " dropped"))

# =============================================================================
# Write output files and generate summary
# =============================================================================

message("\n Writing output files...")

# Write out filtered files
write_tsv(sqanti_df, file.path(dir, paste0(basename, "_classification_filtered.txt")))

save_filtered_fasta(
  sqanti_fasta, 
  kept_ids, 
  file.path(dir, paste0(basename, "_corrected_filtered.fasta"))
)

save_filtered_gtf(
  sqanti_gtf, 
  kept_ids, 
  file.path(dir, paste0(basename, "_corrected_filtered.gtf"))
)

# Filtered hashid and cpm table
all_ids %>% filter(isoform_id %in% sqanti_df$isoform) %>%
  write_tsv(file.path(dir, paste0(basename, "_hashids_with_cpm_filtered.txt")))

# Create dropout reasons table and write dropout files
dropout_reasons_df = map_dfr(names(dropout_tracker), ~ {
  tibble(
    isoform        = dropout_tracker[[.x]],
    dropout_reason = .x
    )
})

dropout_reasons_df %<>%
  left_join(sqanti_df_full)

write_tsv(dropout_reasons_df, file.path(dropout_dir, paste0(basename, "_dropout_transcripts.tsv")))

save_filtered_fasta(
  sqanti_fasta,
  all_dropout_ids,
  file.path(dropout_dir, paste0(basename, "_corrected_dropout.fasta"))
)

save_filtered_gtf(
  sqanti_gtf,
  all_dropout_ids,
  file.path(dropout_dir, paste0(basename, "_corrected_dropout.gtf"))
)

# Print summary
cat("\n=== FILTERING SUMMARY ===\n")
cat(paste0("Sample: ", basename, "\n"))
cat(paste0("Original transcripts: ", length(original_ids), "\n"))
cat(paste0("Final transcripts: ", length(kept_ids), "\n"))
cat(paste0("Retention rate: ", round(100 * length(kept_ids) / length(original_ids), 1), "%\n"))
cat("\nDropout breakdown:\n") 

# Print dropout breakdown
for (reason in names(dropout_tracker)) {
  count = length(dropout_tracker[[reason]])
  pct   = round(100 * count / length(original_ids), 1)
  cat(paste0("  ", reason, ": ", count, " (", pct, "%)\n"))
}