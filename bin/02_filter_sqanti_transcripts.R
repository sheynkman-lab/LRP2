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
  library(optparse)
})

options(scipen = 999)

# =============================================================================
# Get command line arguments and check required files
# =============================================================================

option_list = list(
  make_option(c("--basename"), type = "character", default = NULL,
              help = "Output base name"),
  make_option(c("--classification"), type = "character", default = NULL,
              help = "Path to SQANTI classification file"),
  make_option(c("--sample_gtf"), type = "character", default = NULL,
              help = "Path to SQANTI corrected GTF"),
  make_option(c("--sample_fasta"), type = "character", default = NULL,
              help = "Path to SQANTI corrected FASTA"),
  make_option(c("--mapping_file"), type = "character", default = NULL,
              help = "Path to hashids mapping file from 00_generate_hashids.R"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Output directory for results"),
  make_option(c("--filter_protein_coding"), type = "logical", default = TRUE,
              help = "Filter to protein coding genes only [default: %default]"),
  make_option(c("--filter_internal_priming"), type = "logical", default = TRUE,
              help = "Filter internal priming artifacts [default: %default]"),
  make_option(c("--filter_rts"), type = "logical", default = TRUE,
              help = "Filter template switching artifacts [default: %default]"),
  make_option(c("--percent_polya_threshold"), type = "double", default = 60,
              help = "Percent polyA threshold for internal priming [default: %default]"),
  make_option(c("--transcript_class_keep"), type = "character", default = "FSM,ISM,NIC,NNC",
              help = "Comma-separated transcript classes to keep [default: %default]")
)

opt = parse_args(OptionParser(option_list = option_list))

required_args = c("basename", "classification", "sample_gtf", "sample_fasta",
                  "mapping_file", "output_dir")
missing = required_args[sapply(required_args, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste0("--", missing, collapse = ", "))
}

basename                = opt$basename
classification_file     = opt$classification
sqanti_gtf              = opt$sample_gtf
sqanti_fasta            = opt$sample_fasta
mapping_file            = opt$mapping_file
output_dir              = opt$output_dir
filter_protein_coding   = opt$filter_protein_coding
filter_internal_priming = opt$filter_internal_priming
filter_RTS              = opt$filter_rts
percent_polyA_threshold = opt$percent_polya_threshold
tclass_to_keep          = opt$transcript_class_keep

stopifnot("Classification file not found" = file.exists(classification_file))
stopifnot("Sample GTF not found"          = file.exists(sqanti_gtf))
stopifnot("Sample FASTA not found"        = file.exists(sqanti_fasta))
stopifnot("Mapping file not found"        = file.exists(mapping_file))

# =============================================================================
# Helper functions
# =============================================================================

#' Convert GTF to colored BED12 for genome browser visualization
#'
#' @param gtf_path Path to input GTF file
#' @param output_bed Path to output BED12 file
#' @param color_by Column name for coloring (numeric, 0-1). High values = purple, low = yellow. Default: NULL (black)
gtf_to_bed12 <- function(gtf_path, output_bed, color_by = NULL) {
  
  # Import GTF
  gtf <- import(gtf_path)
  df <- as.data.frame(gtf)
  dt <- as.data.table(df)
  
  # Get transcripts
  tx <- dt[type == "transcript"]
  
  # Determine what to use for blocks
  has_exon <- any(dt$type == "exon")
  has_cds <- any(dt$type == "CDS")
  
  if(has_exon) {
    ex <- dt[type == "exon"]      # Use exons for blocks (shows UTRs)
  } else if(has_cds) {
    ex <- dt[type == "CDS"]       # Use CDS for blocks (CDS only)
  } else {
    stop("GTF must contain either exon or CDS features")
  }
  
  cds <- dt[type == "CDS"]        # Always get CDS for thick regions
  
  # Set keys for fast joins
  setkey(ex, transcript_id)
  if(nrow(cds) > 0) setkey(cds, transcript_id)
  
  # Sort transcripts by gene, then by ratio (descending)
  if(!is.null(color_by)) {
    tx[, sort_val := as.numeric(get(color_by))]
    tx <- tx[order(gene_id, -sort_val)]
  }
  
  # Pre-compute colors
  if(!is.null(color_by)) {
    vals <- as.numeric(tx[[color_by]])
    col_idx <- round((1 - vals) * 99) + 1
    col_idx[is.na(col_idx)] <- 50
    hex_colors <- viridis::viridis(100)[col_idx]
    rgb_matrix <- col2rgb(hex_colors)
    tx$color <- paste(rgb_matrix[1,], rgb_matrix[2,], rgb_matrix[3,], sep=",")
  } else {
    tx$color <- "0,0,0"
  }
  
  # Build BED12 rows
  bed_list <- vector("list", nrow(tx))
  
  for(i in 1:nrow(tx)) {
    t <- tx[i]
    e <- ex[.(t$transcript_id)][order(start)]
    
    if(nrow(e) == 0) next
    
    # Get CDS boundaries for thick regions
    if(nrow(cds) > 0) {
      t_cds <- cds[.(t$transcript_id)]
      if(nrow(t_cds) > 0) {
        thick_start <- min(t_cds$start) - 1
        thick_end <- max(t_cds$end)
      } else {
        # Has exons but no CDS = all thin (UTR only)
        thick_start <- t$start - 1
        thick_end <- t$start - 1
      }
    } else {
      # No CDS at all = all thin
      thick_start <- t$start - 1
      thick_end <- t$start - 1
    }
    
    # Build block strings
    block_sizes <- paste(e$width, collapse=",")
    block_starts <- paste(e$start - t$start, collapse=",")
    
    # Build BED12 line
    bed_list[[i]] <- c(
      as.character(t$seqnames), 
      as.character(t$start-1), 
      as.character(t$end), 
      as.character(t$name), 
      "0", 
      as.character(t$strand),
      as.character(thick_start), 
      as.character(thick_end), 
      as.character(t$color), 
      as.character(nrow(e)),
      block_sizes,
      block_starts
    )
  }
  
  # Remove NULLs and write
  bed <- do.call(rbind, bed_list[!sapply(bed_list, is.null)])
  write.table(bed, output_bed, sep="\t", quote=F, row.names=F, col.names=F)
  
}

# =================================================================================
# Read mapping file and SQANTI classification, calculate CPM
# =================================================================================

# Mapping file from 00_generate_hashids.R (has isoform_id, gene_name, gene_type, etc.)
mapping = read_tsv(mapping_file, show_col_types = FALSE)

# SQANTI classification file- add sample names, calculate cpm on all transcripts before filtering
sqanti_df = read_tsv(classification_file, show_col_types = FALSE)

colnames(sqanti_df) = ifelse(grepl("^FL\\.", colnames(sqanti_df)),
                             paste0(sub("^FL\\.", "", colnames(sqanti_df)), "_counts"),
                             colnames(sqanti_df))

counts_cols = grep("_counts$", colnames(sqanti_df), value = TRUE)
sqanti_df %<>%
  mutate(across(all_of(counts_cols), 
                ~ (.x / sum(.x)) * 1e6,
                .names = "{gsub('_counts', '_cpm', .col)}"))

# Save original isoform IDs before replacement
original_isoform_ids = sqanti_df$isoform

# Replace original isoform IDs with new pipeline isoform_id throughout
sqanti_df_full = sqanti_df %>%
  inner_join(mapping %>% select(isoform_id, original_transcript_id,
                                any_of("gene_type"), any_of("gene_name")),
             by = c("isoform" = "original_transcript_id")) %>%
  mutate(isoform = isoform_id) %>%
  select(-isoform_id)

# Build the hashids + CPM table using new isoform IDs
all_ids = mapping %>%
  inner_join(sqanti_df_full %>% select(isoform,
                                       ends_with("_counts"),
                                       ends_with("_cpm")),
             by = c("isoform_id" = "isoform"))

write_tsv(all_ids, file.path(output_dir, paste0(basename, ".transcriptome.all_hashids_with_cpm.txt")))

# =============================================================================
# Apply filters sequentially
# =============================================================================

cat("\nStarting filtering of SQANTI output:\n")

# Initialize dropout tracking
# Note: sqanti_df_full$isoform now contains the NEW pipeline IDs
all_new_ids = sqanti_df_full$isoform
sqanti_df_full$dropout_reason = "kept"  # Initialize all as kept

# Keep only transcripts of protein coding genes
sqanti_df = sqanti_df_full
if (filter_protein_coding) {
  ids = sqanti_df$isoform
  sqanti_df %<>% 
    filter(gene_type == "protein_coding")
  
  kept_ids    = sqanti_df$isoform
  dropped_ids = setdiff(ids, kept_ids)
  sqanti_df_full$dropout_reason[sqanti_df_full$isoform %in% dropped_ids] = "not_protein_coding"
  #dropout_tracker[["not_protein_coding"]] = dropped_ids
  
  cat("  Protein coding filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts\n")
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
  sqanti_df_full$dropout_reason[sqanti_df_full$isoform %in% dropped_ids] = "internal_priming"
  #dropout_tracker[["internal_priming"]] = dropped_ids
  
  cat("  Internal priming filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts\n")
}

# Filter out template switching, keep RTS_stage = FALSE (no RTS)
if (filter_RTS) {
  ids = sqanti_df$isoform
  
  count_cols = grep("_counts$", colnames(sqanti_df), value = TRUE)
  if (length(count_cols) > 0) {
    sqanti_df$avg_counts = rowMeans(sqanti_df[, count_cols], na.rm = TRUE)
  } else {
    sqanti_df$avg_counts = NA
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
  sqanti_df_full$dropout_reason[sqanti_df_full$isoform %in% dropped_ids] = "template_switching"

  cat("  Template switching filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts\n")
}

# Apply structural category filter
if (toupper(tclass_to_keep) == "ALL") {
  # Keep everything - no filtering
  ids      = sqanti_df$isoform
  kept_ids = sqanti_df$isoform
  cat("No structural category filtering applied\n")
  
} else {
  # track ids before filtering
  ids = sqanti_df$isoform
  
  # Parse comma-separated list
  tclass_list = strsplit(tclass_to_keep, ",")[[1]] %>% 
    trimws()
  
  category_map <- c(
    "FSM" = "full-splice_match",
    "ISM" = "incomplete-splice_match", 
    "NIC" = "novel_in_catalog",
    "NNC" = "novel_not_in_catalog"
  )
  
  allowed_categories = tclass_list %>%
    map_chr(~category_map[.x] %||% .x)
  
  cat("\nFiltering to categories: ", paste(allowed_categories, collapse = ", "), "\n")
  
  sqanti_df %<>% 
    filter(structural_category %in% allowed_categories)
  
  kept_ids = sqanti_df$isoform
}
  
dropped_ids = setdiff(ids, kept_ids)
sqanti_df_full$dropout_reason[sqanti_df_full$isoform %in% dropped_ids] = "structural_category"
cat("  Structural category filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts\n")

# Calculate final results
all_dropout_ids = sqanti_df_full$isoform[sqanti_df_full$dropout_reason != "kept"]

# =============================================================================
# Write output files and generate summary
# =============================================================================

cat("\nWriting output files:\n")

# Build lookup from new isoform_id to original transcript_id
id_lookup = mapping %>% select(isoform_id, original_transcript_id)

#write_tsv(sqanti_df_full, file.path(output_dir, paste0(basename, ".transcriptome.classification_all.txt")))
write_tsv(sqanti_df, file.path(output_dir, paste0(basename, ".transcriptome.filtered_SQANTI_classification.txt")))

# FASTA — filter to kept transcripts and rename headers to new isoform_id
sequences = readDNAStringSet(sqanti_fasta)
kept_original_ids = id_lookup %>% filter(isoform_id %in% kept_ids) %>% pull(original_transcript_id)
filtered_sequences = sequences[names(sequences) %in% kept_original_ids]

# Create name mapping: original_id -> new_id
# Need to swap column order so deframe() creates the correct mapping direction
name_map = id_lookup %>%
  filter(original_transcript_id %in% names(filtered_sequences)) %>%
  select(original_transcript_id, isoform_id) %>%  # Swap order: original first, new second
  deframe()

# Check for any sequences that don't have a mapping
unmapped = setdiff(names(filtered_sequences), names(name_map))
if (length(unmapped) > 0) {
  warning("Found ", length(unmapped), " sequences without isoform_id mapping. Removing them.")
  cat("\nUnmapped transcripts: ", paste(head(unmapped, 10), collapse = ", "), "\n")
  filtered_sequences = filtered_sequences[names(filtered_sequences) %in% names(name_map)]
}

# Apply the name mapping
names(filtered_sequences) = name_map[names(filtered_sequences)]

# Final check for NAs
if (any(is.na(names(filtered_sequences)))) {
  stop("ERROR: Found NA values in sequence names after mapping. This should not happen.")
}

writeXStringSet(filtered_sequences, file.path(output_dir, paste0(basename, ".transcriptome.filtered.fasta")))
cat("  Saved ", length(filtered_sequences), " sequences to ", basename, ".transcriptome.filtered.fasta\n")


# sample gtf
gtf = import(sqanti_gtf) %>%
  as.data.frame()

# Filter to kept transcripts (using original IDs since GTF still has them)
filtered_gtf = gtf %>%
  filter(transcript_id %in% kept_original_ids) %>%
  left_join(id_lookup, by = c("transcript_id" = "original_transcript_id")) %>%
  mutate(transcript_id = isoform_id) %>%
  select(-isoform_id)

# Add avg_ratio for browser visualization
new_attributes = all_ids %>%
  filter(isoform_id %in% kept_ids) %>%
  mutate(avg_cpm = rowMeans(across(ends_with("_cpm")), na.rm = TRUE)) %>%
  group_by(reference_gene_id) %>%
  mutate(
    gene_total_cpm = sum(avg_cpm, na.rm = TRUE),
    avg_ratio = round(avg_cpm / gene_total_cpm, 3)
  ) %>%
  ungroup() %>%
  select(transcript_id = isoform_id, avg_ratio)

filtered_gtf %<>% 
  left_join(new_attributes, by = c("transcript_id")) %>%
  mutate(name = paste0(transcript_id, "|", avg_ratio))

gr_updated = makeGRangesFromDataFrame(filtered_gtf, keep.extra.columns = TRUE)
gtf_output = file.path(output_dir, paste0(basename, ".transcriptome.filtered.gtf"))
export(gr_updated, gtf_output, format = "gtf")

# convert to bed12
gtf_to_bed12(gtf_path = gtf_output,
             output_bed = file.path(output_dir, paste0(basename, ".transcriptome.filtered.bed")),
             color_by = "avg_ratio")

# Filtered hashid and cpm table
all_ids %>% filter(isoform_id %in% kept_ids) %>%
  write_tsv(file.path(output_dir, paste0(basename, ".transcriptome.filtered_hashids_with_cpm.txt")))

# Print summary
cat("\n=== FILTERING SUMMARY ===\n")
cat(paste0("Sample: ", basename, "\n"))
cat(paste0("Original transcripts: ", length(all_new_ids), "\n"))
cat(paste0("Final transcripts: ", length(kept_ids), "\n"))
cat(paste0("Retention rate: ", round(100 * length(kept_ids) / length(all_new_ids), 1), "%\n"))