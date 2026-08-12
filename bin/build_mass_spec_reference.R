#!/usr/bin/env Rscript

#' Build references for mass spec searching
#' 
#' Constructs a combined protein FASTA and reference table for proteogenomic
#' database searching (e.g. Fragpipe). Supports multiple input modes:
#'
#'   1. GENCODE protein FASTA only
#'   2. LRP2-derived/custom protein FASTA (+/- GENCODE)
#'
#' When a count matrix is provided, LRP/custom transcripts are filtered to
#' only those with CPM > 0 in the specified sample. GENCODE entries are
#' always retained in full.
#'
#' Outputs:
#'   - {prefix}.proteomics.reference.fasta  Standardized headers: transcript_id|gene_id|pclass|status|reference_type
#'   - {prefix}.proteomics.reference.tsv    Matching metadata table
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
# --custom_fasta = custom + GENCODE (concatenated)
# --custom_fasta --no_gencode = custom only


option_list <- list(
  make_option(c("--custom_fasta"), type = "character", default = NULL,
              help = "Path to user-provided protein FASTA (LRP pipeline output or custom)"),
  make_option(c("--gencode_fasta"), type = "character", default = NULL,
              help = "Path to GENCODE protein FASTA (will download if not provided)"),
  make_option(c("--no_gencode"), action = "store_true", default = FALSE,
              help = "Do not include GENCODE in the reference"),
  make_option(c("--sample_name"), type = "character", default = NULL,
              help = "Sample name for output file naming [required]"),
  make_option(c("--counts"), type = "character", default = NULL,
              help = "Path to count matrix TSV (used when --custom_fasta is provided)"),
  make_option(c("--outdir"), type = "character", default = ".",
              help = "Output directory [default: .]")
)

opt = parse_args(OptionParser(option_list = option_list))

# =============================================================================
# Validate inputs
# =============================================================================

if (is.null(opt$sample_name)) stop("--sample_name is required")

if (is.null(opt$custom_fasta) && opt$no_gencode) {
  stop("--no_gencode cannot be used without --custom_fasta. There would be no reference.")
}

has_fasta = !is.null(opt$custom_fasta)

if (has_fasta && is.null(opt$counts)) {
  warning("--counts not provided. Sample specific filtering will be skipped.")
}

if (!is.null(opt$counts) && !has_fasta) {
  warning("--counts provided but no --custom_fasta. Count matrix will be ignored.")
}

if (!opt$no_gencode && (is.null(opt$gencode_fasta) || !file.exists(opt$gencode_fasta))) {
  stop("--gencode_fasta is required unless --no_gencode is set.")
}

# =============================================================================
# Output prefixes (set later after we know what's actually included)
# =============================================================================
# Output prefix will be set after building combined_ref so we know what was actually included


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
gencode_label = NULL

if (!opt$no_gencode) {
  
  cat("Using GENCODE protein FASTA:", opt$gencode_fasta, "\n")
  
  # Label for output naming, taken from the FASTA actually used
  gencode_label = str_extract(basename(opt$gencode_fasta), "vM?\\d+")
  if (is.na(gencode_label)) gencode_label = "gencode"
  
  # Header: ENSP|ENST|ENSG|OTTHUMG|OTTHUMT|TXNAME|GENE|LENGTH
  gencode_raw = fasta_to_tibble(opt$gencode_fasta)
  gencode_ref = gencode_raw %>%
    separate_wider_delim(header,
                         delim = "|",
                         names = c(NA, "transcript_id", "gene_id", NA, NA, NA, NA, NA),
                         cols_remove = FALSE) %>%
    mutate(pclass = NA_character_, status = NA_character_, reference_type = "gencode") %>%
    select(header, transcript_id, gene_id, pclass, status, reference_type, sequence)
  
  combined_ref = bind_rows(combined_ref, gencode_ref)
}

# This chunk auto-detects custom versus LRP2 FASTA and creates the header accordingly
if (!is.null(opt$custom_fasta)) {
  custom_raw = fasta_to_tibble(opt$custom_fasta)
  
  # Validate: headers must contain at least one | separating transcript_id and gene_id
  bad_headers = custom_raw %>% filter(!str_detect(header, "\\|"))
  if (nrow(bad_headers) > 0) {
    stop("Input FASTA headers must contain transcript_id and gene_id separated by '|'. ",
         nrow(bad_headers), " headers lack a '|' delimiter. First example: ", bad_headers$header[1])
  }
  
  # Validate: transcript_ids must be unique
  transcript_ids = str_extract(custom_raw$header, "^[^|]+")
  dupes = transcript_ids[duplicated(transcript_ids)]
  if (length(dupes) > 0) {
    stop("Input FASTA contains ", length(dupes), " duplicate transcript IDs. Check that header is transcript_id | gene_id")
  }
  
  parsed = custom_raw %>%
    separate_wider_delim(header, delim = "|",
                         names = c("transcript_id", "gene_id", "gene_name", "pclass", "status"),
                         too_few = "align_start", too_many = "drop",
                         cols_remove = FALSE)
  
  # LRP2 IDs are gene_name::isoform_id; FPM is an LRP-specific pclass
  is_lrp = all(str_detect(parsed$transcript_id, fixed("::"))) &
    any(parsed$pclass == "FPM", na.rm = TRUE)
  
  if (is_lrp) {
    cat("FASTA source: lrp (transcript IDs contain '::' and FPM found in pclass)\n")
    cat("  Retaining pclass and status from header\n")
  } else {
    cat("FASTA source: custom\n")
    cat("  Only transcript_id and gene_id will be used; pclass and status set to NA\n")
    parsed %<>% mutate(pclass = NA_character_, status = NA_character_)
  }
  
  custom_ref = parsed %>%
    mutate(reference_type = if_else(is_lrp, "lrp", "custom")) %>%
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

  # If no matching CPM column, remove LRP/custom entries and use GENCODE only
  if (length(cpm_cols) == 0) {
    warning("No column found containing both '", opt$sample_name, "' and 'cpm' in count matrix. ",
            "This sample has no matching RNA data. Only the GENCODE reference will be used.")

    # Remove all LRP/custom entries, keep only GENCODE transcripts
    n_before = sum(combined_ref$reference_type != "gencode")
    combined_ref = combined_ref %>% filter(reference_type == "gencode")
    cat("No matching RNA sample found for '", opt$sample_name, "'\n")
    cat("  Removed ", n_before, " LRP/custom entries\n")
    cat("  Using GENCODE reference only (", nrow(combined_ref), " entries)\n")

  } else {
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
}

# =============================================================================
# Write out reference fasta and matching table
# =============================================================================

# Determine output prefix based on what's actually in the combined reference
prefix_parts = c(opt$sample_name)
has_lrp_in_ref = any(combined_ref$reference_type == "lrp")
has_custom_in_ref = any(combined_ref$reference_type == "custom")

if (has_lrp_in_ref)    prefix_parts = c(prefix_parts, "lrp")
if (has_custom_in_ref) prefix_parts = c(prefix_parts, "custom")
if (!opt$no_gencode)   prefix_parts = c(prefix_parts, gencode_label)

out_prefix = paste(prefix_parts, collapse = ".")

cat("\nOutput prefix:", out_prefix, "\n")

# Combined reference table (without sequences)
write_tsv(select(combined_ref, -sequence), file.path(opt$outdir, paste0(out_prefix, ".proteomics.reference.tsv")))

# Combined FASTA with standardized headers: transcript_id|gene_id|pclass|status|reference_type
combined_ref %<>%
  mutate(new_header = paste(transcript_id, gene_id, pclass, status, reference_type, sep = "|"))

fasta_out        = AAStringSet(combined_ref$sequence)
names(fasta_out) = combined_ref$new_header

writeXStringSet(fasta_out, file.path(opt$outdir, paste0(out_prefix, ".proteomics.reference.fasta")))

cat("\n--- Summary ---\n")
cat("Combined reference:", nrow(combined_ref), "total entries\n")
combined_ref %>%
  count(reference_type) %>%
  mutate(label = paste0("  ", reference_type, ": ", n)) %>%
  pull(label) %>%
  walk(cat, "\n")

cat("Done.\n")