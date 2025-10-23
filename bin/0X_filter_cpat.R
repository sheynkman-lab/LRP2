#!/usr/bin/env Rscript

#' CPAT Coding Potential Filter
#' 
#' Filters transcripts based on CPAT coding probability scores.
#' Separates likely protein-coding transcripts from non-coding RNAs.
#'
#' Outputs both filtered files and dropout files.

# =============================================================================
# Load required libraries
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(magrittr)
})

# =============================================================================
# Get environment variables
# =============================================================================

basename <- Sys.getenv("OUTPUT_BASE_NAME")
sqanti_dir <- file.path(Sys.getenv("OUTPUT_DIR"), "sqanti") 
cpat_dir <- file.path(Sys.getenv("OUTPUT_DIR"), "orf_calling") 
coding_threshold <- as.numeric(Sys.getenv("CPAT_CODING_THRESHOLD", "0.364"))

# =============================================================================
# Check required files
# =============================================================================

input_fasta_file <- file.path(sqanti_dir, paste0(basename, "_corrected_filtered.fasta"))
stopifnot("Input FASTA file not found" = file.exists(input_fasta_file))

cpat_results_file <- file.path(cpat_dir, paste0(basename, ".ORF_prob.tsv"))
stopifnot("CPAT results file not found" = file.exists(cpat_results_file))

# Create output directories
dropout_dir <- file.path(cpat_dir, "dropout")
dir.create(dropout_dir, recursive = TRUE, showWarnings = FALSE)

message("Starting CPAT-based coding potential filtering...")
message(paste0("Using coding threshold: ", coding_threshold))

# =============================================================================
# LOAD AND FILTER CPAT RESULTS
# =============================================================================

# Load CPAT results
cpat_df <- read_tsv(cpat_results_file, show_col_types = FALSE)

# Apply coding probability filter
coding_df <- cpat_df %>% filter(Coding_prob >= coding_threshold)
dropout_df <- cpat_df %>% filter(Coding_prob < coding_threshold)

# Extract transcript IDs
coding_ids <- coding_df$ID
dropout_ids <- dropout_df$ID

# =============================================================================
# SAVE FILTERED TSV FILES
# =============================================================================

# Save filtered CPAT results
write_tsv(coding_df, file.path(cpat_dir, paste0(basename, "_cpat_filtered.tsv")))
write_tsv(dropout_df, file.path(dropout_dir, paste0(basename, "_cpat_dropout.tsv")))

# =============================================================================
# FILTER FASTA FILES
# =============================================================================

message("Filtering FASTA sequences...")

# Read input FASTA
sequences <- readDNAStringSet(input_fasta_file)

# Filter sequences for coding transcripts
coding_mask <- names(sequences) %in% coding_ids
coding_sequences <- sequences[coding_mask]

# Filter sequences for dropout transcripts  
dropout_mask <- names(sequences) %in% dropout_ids
dropout_sequences <- sequences[dropout_mask]

# Write filtered FASTA files
writeXStringSet(coding_sequences, 
                file.path(cpat_dir, paste0(basename, "_cpat_filtered.fasta")))

writeXStringSet(dropout_sequences,
                file.path(dropout_dir, paste0(basename, "_cpat_dropout.fasta")))

# =============================================================================
# SUMMARY REPORT
# =============================================================================

message("CPAT filtering completed successfully!")

cat("\n=== CPAT FILTERING SUMMARY ===\n")
cat(paste0("Sample: ", basename, "\n"))
cat(paste0("Coding threshold: ", coding_threshold, "\n"))
cat(paste0("Total ORFs: ", nrow(cpat_df), "\n"))
cat(paste0("Likely coding: ", length(coding_ids), "\n"))
cat(paste0("Non-coding: ", length(dropout_ids), "\n"))
cat(paste0("Coding retention rate: ", round(100 * length(coding_ids) / nrow(cpat_df), 1), "%\n"))