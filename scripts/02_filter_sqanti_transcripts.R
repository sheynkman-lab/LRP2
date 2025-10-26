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
  all    = c("antisense", "novel_not_in_catalog", "novel_in_catalog",
             "incomplete-splice_match", "full-splice_match", "genic",
             "intergenic", "fusion", "genic_intron")
)

# =============================================================================
# Get environment variables for filtering parameters
# =============================================================================

filter_protein_coding   <- as.logical(Sys.getenv("PROTEIN_CODING_FILTER", "TRUE"))
filter_internal_priming <- as.logical(Sys.getenv("INTERNAL_PRIMING_FILTER", "TRUE"))
filter_RTS              <- as.logical(Sys.getenv("TEMPLATE_SWITCHING_FILTER", "TRUE"))
percent_polyA_threshold <- as.numeric(Sys.getenv("PERCENT_POLYA_THRESHOLD", "60"))
structural_level        <- Sys.getenv("STRUCTURE_FILTER", "strict")

basename    <- Sys.getenv("OUTPUT_BASE_NAME")
dir         <- file.path(Sys.getenv("OUTPUT_DIR"), "sqanti") 
gencode_gtf <- Sys.getenv("GENCODE_GTF_FILE")

# =============================================================================
# Check required files
# =============================================================================

classification_file = paste0(dir, "/", basename, "_classification.txt")
sqanti_gtf          = paste0(dir, "/", basename, "_corrected.gtf")
sqanti_fasta        = paste0(dir, "/", basename, "_corrected.fasta")

stopifnot("File not found" = file.exists(classification_file))
stopifnot("File not found" = file.exists(sqanti_gtf))
stopifnot("File not found" = file.exists(sqanti_fasta))

# Create additional output directory
dropout_dir = file.path(dir, "dropout")
dir.create(dropout_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Helper functions
# =============================================================================

#' Save filtered sequences to FASTA
#' @param fasta_file Input FASTA file
#' @param keep_ids IDs to keep
#' @param output_file Output FASTA file
save_filtered_fasta <- function(fasta_file, keep_ids, output_file) {

  sequences          <- readDNAStringSet(fasta_file)
  keep_mask          <- names(sequences) %in% keep_ids
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

  gtf          <- import(gtf_file)
  filtered_gtf <- gtf[gtf$transcript_id %in% keep_ids]
  
  # Write filtered GTF
  export(filtered_gtf, output_file)
  message(paste0("Saved ", length(unique(filtered_gtf$transcript_id)), " transcripts to ", basename(output_file)))
}

# =============================================================================
# Load and merge inputs
# =============================================================================

message("\n Starting filtering of SQANTI output...")

# Gencode gtf used to select protein coding genes and add column of readable gene name
gencode    = import(gencode_gtf, format = "gtf")
gencode_df = as.data.table(gencode)

gencode_df %<>% 
  filter(type == "transcript") %>% 
  select(associated_gene = gene_id, gene_type, gene_name, associated_transcript = transcript_id, transcript_name, protein_id)

gencode_gene = gencode_df %>% 
  select(associated_gene, gene_type, gene_name) %>% 
  distinct()

# SQANTI classification file- keep only isoforms with ENS gene IDs
sqanti_df_full = read_tsv(classification_file, show_col_types = FALSE) %>%
  filter(!is.na(associated_gene), str_starts(associated_gene, "ENS"))

sqanti_df_full %<>% 
  left_join(gencode_gene) %>% 
  left_join(gencode_df)


# =============================================================================
# Apply filters sequentially
# =============================================================================

# Initialize dropout tracking
original_ids    = sqanti_df_full$isoform
dropout_tracker = list()

# Keep only transcripts of protein coding genes
sqanti_df = sqanti_df_full
if (filter_protein_coding) {
  ids = sqanti_df$isoform
  sqanti_df %<>% 
    filter(gene_type == "protein_coding")
  
  kept_ids    = sqanti_df$isoform
  dropped_ids = setdiff(ids, kept_ids)
  dropout_tracker[["not_protein_coding"]] = dropped_ids
  
  message("Protein coding filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts")
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
  if (gene_type_pos - ratio_tss_pos > 1) {
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

# Create dropout reasons table
dropout_reasons_df = map_dfr(names(dropout_tracker), ~ {
  tibble(
    isoform        = dropout_tracker[[.x]], 
    dropout_reason = .x
    )
})

dropout_reasons_df %<>% 
  left_join(sqanti_df_full)

# Write out dropout files
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