#!/usr/bin/env Rscript

#' Build references for mass spec searching
#' 
#' Constructs a combined protein FASTA and reference table for proteogenomic
#' database searching (e.g. Fragpipe). Supports three input modes:
#'
#'   1. GENCODE only
#'   2. LRP pipeline FASTA (+/- GENCODE)
#'   3. Custom user FASTA (+/- GENCODE)
#'
#' When a count matrix is provided, LRP/custom transcripts are filtered to
#' only those with CPM > 0 in the specified sample. GENCODE entries are
#' always retained in full.
#'
#' Outputs:
#'   - {prefix}_reference.fasta  Standardized headers: transcript_id|gene_id|pclass|status|reference_type
#'   - {prefix}_reference.tsv    Matching metadata table
#'
#' Usage:
#'   Rscript build_mass_spec_reference.R --sample_name sample --lrp_fasta lrp_orfs.fasta --gencode_version 47 --counts counts.tsv --outdir results/
#'   


# =============================================================================
# Load libraries
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(Biostrings)
  library(optparse)
})

options(scipen = 999)


# =============================================================================
# Command Line Arguments
# =============================================================================

# The logic is

# No flags = downloads and uses GENCODE only
# --lrp_fasta = LRP + GENCODE (concatenated)
# --custom_fasta = custom + GENCODE (concatenated)
# --lrp_fasta --no_gencode = LRP only
# --custom_fasta --no_gencode = custom only


option_list <- list(
  make_option(c("--lrp_fasta"), type = "character", default = NULL,
              help = "Path to lrp pipeline protein FASTA"),
  make_option(c("--custom_fasta"), type = "character", default = NULL,
              help = "Path to user provided protein FASTA."),
  make_option(c("--gencode_version"), type = "character", default = "46",
              help = "GENCODE version to download [default: 46]"),
  make_option(c("--no_gencode"), action = "store_true", default = FALSE,
              help = "Do not include GENCODE in the reference"),
  make_option(c("--sample_name"), type = "character", default = NULL,
              help = "Sample name for output file naming [required]"),
  make_option(c("--counts"), type = "character", default = NULL,
              help = "Path to count matrix TSV (required when --lrp_fasta or --custom_fasta is provided)"),
  make_option(c("--outdir"), type = "character", default = ".",
              help = "Output directory [default: .]")
)

opt = parse_args(OptionParser(option_list = option_list))

# =============================================================================
# Validate inputs
# =============================================================================

#if (is.null(opt$counts)) stop("--counts is required")
if (is.null(opt$sample_name)) stop("--sample_name is required")
if (!is.null(opt$lrp_fasta) && !is.null(opt$custom_fasta)) {
  stop("--lrp_fasta and --custom_fasta are mutually exclusive. Provide one or neither.")
}
if (is.null(opt$lrp_fasta) && is.null(opt$custom_fasta) && opt$no_gencode) {
  stop("--no_gencode cannot be used without --lrp_fasta or --custom_fasta. There would be no reference.")
}

has_fasta = !is.null(opt$lrp_fasta) || !is.null(opt$custom_fasta)
if (has_fasta && is.null(opt$counts)) {
  warning("--counts not provided. Sample specific filtering will be skipped.")
}
if (!is.null(opt$counts) && !has_fasta) {
  warning("--counts provided but no --lrp_fasta or --custom_fasta. Count matrix will be ignored.")
}

# =============================================================================
# Output prefixes
# =============================================================================

prefix_parts = c()
if (!is.null(opt$lrp_fasta))    prefix_parts = c(prefix_parts, paste0(opt$sample_name, ".lrp"))
if (!is.null(opt$custom_fasta)) prefix_parts = c(prefix_parts, paste0(opt$sample_name, ".custom"))
if (!opt$no_gencode)            prefix_parts = c(prefix_parts, paste0("gencode.v", opt$gencode_version))
out_prefix = paste(prefix_parts, collapse = "_")


# =============================================================================
# Helper functions
# =============================================================================

fasta_to_tibble = function(path) {
  seqs = readAAStringSet(path)
  tibble(header = names(seqs), sequence = as.character(seqs))
}

# =============================================================================
# Build reference
# =============================================================================

combined_ref = tibble()

# GENCODE (included by default unless --no_gencode)
if (!opt$no_gencode) {
  
  gencode_url  = paste0("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_", opt$gencode_version, "/gencode.v", opt$gencode_version, ".pc_translations.fa.gz")
  gencode_path = file.path(opt$outdir, basename(gencode_url))
  
  if (!file.exists(gencode_path)) {
    cat("Downloading GENCODE v", opt$gencode_version, " protein FASTA...\n")
    download.file(gencode_url, gencode_path, mode = "wb")
  }
  
  # Header: ENSP|ENST|ENSG|OTTHUMG|OTTHUMT|TXNAME|GENE|LENGTH
  gencode_raw = fasta_to_tibble(gencode_path)
  gencode_ref = gencode_raw %>%
    separate_wider_delim(header, 
                         delim = "|",
                         names = c(NA, "transcript_id", "gene_id", NA, NA, NA, NA, NA), 
                         cols_remove = FALSE) %>%
    mutate(pclass = "FSM", status = "gencode", reference_type = "gencode") %>%
    select(header, transcript_id, gene_id, pclass, status, reference_type, sequence)
  
  combined_ref = bind_rows(combined_ref, gencode_ref)
}

# LRP pipeline FASTA
if (!is.null(opt$lrp_fasta)) {
  
  # Header: transcript_id|gene_id|gene_name|pclass|status
  lrp_raw = fasta_to_tibble(opt$lrp_fasta)
  lrp_ref = lrp_raw %>%
    separate_wider_delim(header, 
                         delim = "|",
                         names = c("transcript_id", "gene_id", NA, "pclass", "status"),
                         cols_remove = FALSE) %>%
    mutate(reference_type = "lrp") %>%
    select(header, transcript_id, gene_id, pclass, status, reference_type, sequence)
  
  combined_ref = bind_rows(combined_ref, lrp_ref)
  
}

# Custom user FASTA
if (!is.null(opt$custom_fasta)) {
  custom_raw = fasta_to_tibble(opt$custom_fasta)
  
  # Validate: headers must contain at least one | separating transcript_id and gene_id
  bad_headers = custom_raw %>% filter(!str_detect(header, "\\|"))
  if (nrow(bad_headers) > 0) {
    stop("Custom FASTA headers must contain transcript_id and gene_id separated by '|'. ",
         nrow(bad_headers), " headers lack a '|' delimiter. First example: ", bad_headers$header[1])
  }
  
  # Validate: transcript_ids must be unique
  transcript_ids = str_extract(custom_raw$header, "^[^|]+")
  dupes = transcript_ids[duplicated(transcript_ids)]
  if (length(dupes) > 0) {
    stop("Custom FASTA contains ", length(dupes), " duplicate transcript IDs. Check that header is transcript_id | gene_id")
  }
  
  custom_ref = custom_raw %>%
    separate_wider_delim(header, delim = "|",
                         names = c("transcript_id", "gene_id"),
                         too_many = "drop", cols_remove = FALSE) %>%
    mutate(
      pclass         = NA_character_,
      status         = "custom",
      reference_type = "custom"
    ) %>%
    select(header, transcript_id, gene_id, pclass, status, reference_type, sequence)
  
  combined_ref = bind_rows(combined_ref, custom_ref)
  
}

# =============================================================================
# Count matrix - per sample_name in lrp or custom
# =============================================================================

if (has_fasta && !is.null(opt$counts)) {
  counts = read_tsv(opt$counts, show_col_types = FALSE)

  # Select transcript_id (column 1) and CPM column matching sample name
  cpm_cols = names(counts)[str_detect(names(counts), fixed(opt$sample_name)) &
                             str_detect(names(counts), fixed("cpm"))]
  if (length(cpm_cols) == 0) {
    stop("No column found containing both '", opt$sample_name, "' and 'cpm' in count matrix.")
  }
  if (length(cpm_cols) > 1) {
    warning("Multiple CPM columns matched: ", paste(cpm_cols, collapse = ", "), ". Using first: ", cpm_cols[1])
  }
  
  sample_counts = counts %>%
    select(transcript_id = 1, sample_cpm = all_of(cpm_cols[1])) %>%
    filter(sample_cpm > 0)
  
  # Check overlap with FASTA transcript IDs
  fasta_ids = combined_ref %>% filter(reference_type != "gencode") %>% pull(transcript_id)
  
  n_overlap = sum(fasta_ids %in% sample_counts$transcript_id)
  if (n_overlap == 0) {
    stop("No transcript IDs in FASTA match the count matrix. Check that IDs are consistent between files.")
  }
  
  # Filter combined_ref: keep gencode entries + only LRP/custom transcripts with CPM > 0
  n_before = sum(combined_ref$reference_type != "gencode")
  combined_ref %<>%
    filter(reference_type == "gencode" | transcript_id %in% sample_counts$transcript_id) %>%
    left_join(select(sample_counts, transcript_id, sample_cpm), by = "transcript_id")
  n_after = sum(combined_ref$reference_type != "gencode")
  
  cat("Starting LRP/custom FASTA entries:", n_before, "\n")
  cat("  Removed (not expressed in", opt$sample_name, "):", n_before - n_after, "\n")
  cat("  Kept:", n_after, "\n")
  
}

# =============================================================================
# Write out reference fasta and matching table
# =============================================================================

# Combined reference table (without sequences)
write_tsv(select(combined_ref, -sequence), file.path(opt$outdir, paste0(out_prefix, "_reference.tsv")))

# Combined FASTA with standardized headers: transcript_id|gene_id|pclass|status|reference_type
combined_ref %<>%
  mutate(new_header = paste(transcript_id, gene_id, pclass, status, reference_type, sep = "|"))

fasta_out        = AAStringSet(combined_ref$sequence)
names(fasta_out) = combined_ref$new_header

writeXStringSet(fasta_out, file.path(opt$outdir, paste0(out_prefix, "_reference.fasta")))

cat("\n--- Summary ---\n")
cat("Combined reference:", nrow(combined_ref), "total entries\n")
combined_ref %>%
  count(reference_type) %>%
  mutate(label = paste0("  ", reference_type, ": ", n)) %>%
  pull(label) %>%
  walk(cat, "\n")

cat("Done.\n")