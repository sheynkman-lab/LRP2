#!/usr/bin/env Rscript

#' SQANTI3 Transcript Filtering Script
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
  library(tidyverse)
  library(rtracklayer)
  library(Biostrings) 
  library(data.table)
  library(magrittr)
})

# Suppress dplyr informational messages
options(dplyr.inform = FALSE)

# =============================================================================
# Define SQANTI structural categories
# =============================================================================

STRUCTURAL_CATEGORIES <- list(
  strict = c("novel_not_in_catalog", "novel_in_catalog", 
             "incomplete-splice_match", "full-splice_match"),
  all = c("antisense", "novel_not_in_catalog", "novel_in_catalog",
          "incomplete-splice_match", "full-splice_match", "genic",
          "intergenic", "fusion", "genic_intron")
)

# =============================================================================
# Get environment variables for filtering parameters
# =============================================================================

filter_protein_coding <- as.logical(Sys.getenv("PROTEIN_CODING_FILTER", "TRUE"))
filter_internal_priming <- as.logical(Sys.getenv("INTERNAL_PRIMING_FILTER", "TRUE"))
filter_RTS <- as.logical(Sys.getenv("TEMPLATE_SWITCHING_FILTER", "TRUE"))
percent_polyA_threshold <- as.numeric(Sys.getenv("PERCENT_POLYA_THRESHOLD", "95"))
structural_level <- Sys.getenv("STRUCTURE_FILTER", "strict")

basename <- Sys.getenv("OUTPUT_BASE_NAME")
dir <- file.path(Sys.getenv("OUTPUT_DIR"), "sqanti") 
gencode_gtf <- Sys.getenv("GENCODE_GTF_FILE")

# =============================================================================
# Check required files
# =============================================================================

classification_file <- paste0(dir, "/", basename, "_classification.txt")
stopifnot("File not found" = file.exists(classification_file))

sqanti_gtf <- paste0(dir, "/", basename, "_corrected.gtf")
stopifnot("File not found" = file.exists(sqanti_gtf))

sqanti_fasta <- paste0(dir, "/", basename, "_corrected.fasta")
stopifnot("File not found" = file.exists(sqanti_fasta))

# Create additional output directory
dropout_dir <- file.path(dir, "dropout")
dir.create(dropout_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Helper functions
# =============================================================================

#' Apply filter and track dropouts
#' @param df Data frame to filter
#' @param condition Logical condition for filtering
#' @param reason Reason for dropout
#' @param dropout_tracker List to store dropout reasons
#' @return Filtered data frame
apply_filter <- function(df, condition, reason, dropout_tracker) {
  pre_ids <- df$isoform
  df_filtered <- df[condition, ]
  post_ids <- df_filtered$isoform
  
  # Track dropped transcripts
  dropped <- setdiff(pre_ids, post_ids)
  dropout_tracker[[reason]] <- dropped
  
  return(list(
    data = df_filtered,
    tracker = dropout_tracker
  ))
}

#' Save filtered sequences to FASTA
#' @param fasta_file Input FASTA file
#' @param keep_ids IDs to keep
#' @param output_file Output FASTA file
save_filtered_fasta <- function(fasta_file, keep_ids, output_file) {
  # Read FASTA
  sequences <- readDNAStringSet(fasta_file)
  
  # Filter sequences
  keep_mask <- names(sequences) %in% keep_ids
  filtered_sequences <- sequences[keep_mask]
  
  # Write filtered FASTA
  writeXStringSet(filtered_sequences, output_file)
  message(paste0("Saved ", length(filtered_sequences), " sequences to ", basename(output_file)))
}

#' Save filtered GTF
#' @param gtf_file Input GTF file
#' @param keep_ids Transcript IDs to keep
#' @param output_file Output GTF file
save_filtered_gtf <- function(gtf_file, keep_ids, output_file) {
  # Read GTF
  gtf <- import(gtf_file)
  
  # Filter by transcript_id
  filtered_gtf <- gtf[gtf$transcript_id %in% keep_ids]
  
  # Write filtered GTF
  export(filtered_gtf, output_file)
  message(paste0("Saved ", length(unique(filtered_gtf$transcript_id)), " transcripts to ", basename(output_file)))
}

# =============================================================================
# Main code
# =============================================================================

message("\n Loading SQANTI3 classification data...")

# Load and process gencode gtf file
gencode <- import(gencode_gtf, format = "gtf")
gencode_df <- as.data.table(gencode)

gencode_df %<>% filter(type == "transcript") %>% select(gene_id, gene_type, gene_name, transcript_id, transcript_type, transcript_name, protein_id)
gencode_gene = gencode_df %>% select(associated_gene = gene_id, gene_type, gene_name) %>% distinct()

# Load and preprocess classification data, merge with gencode annotation info
sqanti_df_full <- read_tsv(classification_file, show_col_types = FALSE) %>%
  # Remove rows without associated genes and keep only ENS gene IDs
  filter(!is.na(associated_gene), 
         str_starts(associated_gene, "ENS"))

sqanti_df_full %<>% left_join(gencode_gene)

original_ids <- sqanti_df_full$isoform
message(paste0("Starting with ", length(original_ids), " transcripts"))

# Initialize dropout tracking
dropout_tracker <- list()

# Apply filters sequentially
sqanti_df = sqanti_df_full
if (filter_protein_coding) {
  message("Applying protein-coding gene filter...")
  result <- apply_filter(
    sqanti_df,
    sqanti_df$gene_type == "protein_coding",
    "not_protein_coding",
    dropout_tracker
  )
  
  sqanti_df <- result$data
  dropout_tracker <- result$tracker
}

if (filter_internal_priming) {
  message("Applying internal priming filter...")
  result <- apply_filter(
    sqanti_df,
    sqanti_df$perc_A_downstream_TTS <= percent_polyA_threshold,
    "internal_priming",
    dropout_tracker
  )
  
  sqanti_df <- result$data
  dropout_tracker <- result$tracker
}

if (filter_RTS) {
  message("Applying template switching filter...")
  result <- apply_filter(
    sqanti_df,
    !sqanti_df$RTS_stage,  # Keep only FALSE (no RTS)
    "Template_switching_artifact",
    dropout_tracker
  )
  
  sqanti_df <- result$data
  dropout_tracker <- result$tracker
}

# Apply structural category filter
message("Applying structural category filter...")
allowed_categories <- STRUCTURAL_CATEGORIES[[structural_level]]
result <- apply_filter(
  sqanti_df,
  sqanti_df$structural_category %in% allowed_categories,
  "structural_category_filtered",
  dropout_tracker
)

sqanti_df <- result$data
dropout_tracker <- result$tracker

# Calculate final results
kept_ids <- sqanti_df$isoform
all_dropout_ids <- setdiff(original_ids, kept_ids)

message(paste0("Filtering is finished: ", length(kept_ids), " kept, ", length(all_dropout_ids), " dropped"))

# =============================================================================
# Write output files and generate summary
# =============================================================================

message("\n Writing output files...")
write_tsv(sqanti_df, file.path(dir, paste0(basename, "_classification_filtered.txt")))

save_filtered_fasta(sqanti_fasta, kept_ids, file.path(dir, paste0(basename, "_corrected_filtered.fasta")))

save_filtered_gtf(sqanti_gtf, kept_ids, file.path(dir, paste0(basename, "_corrected_filtered.gtf")))

# Save dropout files
dropout_reasons_df <- map_dfr(names(dropout_tracker), ~ {
  tibble(isoform = dropout_tracker[[.x]], dropout_reason = .x)
})

dropout_reasons_df %<>% left_join(sqanti_df_full)

write_tsv(dropout_reasons_df, file.path(dropout_dir, paste0(basename, "_dropout_transcripts.tsv")))
save_filtered_fasta(sqanti_fasta, all_dropout_ids, file.path(dropout_dir, paste0(basename, "_corrected_dropout.fasta")))
save_filtered_gtf(sqanti_gtf, all_dropout_ids, file.path(dropout_dir, paste0(basename, "_corrected_dropout.gtf")))

# Print summary
cat("\n=== FILTERING SUMMARY ===\n")
cat(paste0("Sample: ", basename, "\n"))
cat(paste0("Original transcripts: ", length(original_ids), "\n"))
cat(paste0("Final transcripts: ", length(kept_ids), "\n"))
cat(paste0("Retention rate: ", round(100 * length(kept_ids) / length(original_ids), 1), "%\n"))
cat("\nDropout breakdown:\n") 

dropout_summary <- dropout_reasons_df %>%
  count(dropout_reason, sort = TRUE)

for (i in seq_len(nrow(dropout_summary))) {
  cat(paste0("  ", dropout_summary$dropout_reason[i], ": ", dropout_summary$n[i], "\n"))
}
