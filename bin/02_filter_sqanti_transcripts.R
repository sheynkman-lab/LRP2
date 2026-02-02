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

options(scipen = 999)

# =============================================================================
# Get environment variables for filtering parameters
# =============================================================================

dir                     <- file.path(Sys.getenv("OUTPUT_DIR"), "sqanti_transcript")
basename                <- Sys.getenv("OUTPUT_BASE_NAME")
gencode_gtf             <- Sys.getenv("GENCODE_GTF_FILE")
metadata_path           <- Sys.getenv("SAMPLE_METADATA")
filter_protein_coding   <- as.logical(Sys.getenv("PROTEIN_CODING_FILTER", "TRUE"))
filter_internal_priming <- as.logical(Sys.getenv("INTERNAL_PRIMING_FILTER", "TRUE"))
filter_RTS              <- as.logical(Sys.getenv("TEMPLATE_SWITCHING_FILTER", "TRUE"))
percent_polyA_threshold <- as.numeric(Sys.getenv("PERCENT_POLYA_THRESHOLD", "60"))
tclass_to_keep          <- Sys.getenv("TRANSCRIPT_CLASS_KEEP", "FSM,ISM,NIC,NNC")
#structural_level       <- Sys.getenv("STRUCTURE_FILTER", "strict")

# =============================================================================
# Check required files
# =============================================================================

classification_file = paste0(dir, "/", basename, "_classification.txt")
sqanti_gtf          = paste0(dir, "/", basename, "_corrected.gtf")
sqanti_fasta        = paste0(dir, "/", basename, "_corrected.fasta")
mapping_file        = paste0(dir, "/", basename, "_hashids_mapping.txt")

stopifnot("File not found" = file.exists(classification_file))
stopifnot("File not found" = file.exists(sqanti_gtf))
stopifnot("File not found" = file.exists(sqanti_fasta))
stopifnot("File not found" = file.exists(mapping_file))

# Create additional output directory
#dropout_dir = file.path(dir, "dropout")
#dir.create(dropout_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Helper functions
# =============================================================================

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
# Combine sqanti associated gene_ids with reference ids and hash ids, include counts
# =================================================================================

# reference gtf used to select protein coding genes and add column of readable gene name
gencode    = import(gencode_gtf, format = "gtf")
gencode_df = as.data.frame(gencode)

gencode_df %<>% 
  filter(type == "transcript") %>% 
  select(associated_gene = gene_id, 
         gene_type = any_of(c("gene_type", "gene_biotype")), 
         gene_name = any_of(c("gene_name", "gene")),
         associated_transcript = transcript_id, 
         any_of("transcript_name"), 
         any_of("protein_id"))

gencode_gene = gencode_df %>% 
  select(associated_gene, any_of("gene_type"), any_of("gene_name")) %>% 
  distinct()

# Metatable to add specific sample ids
metadata      = fread(metadata_path)
sample_lookup = setNames(metadata$sample_name, metadata$bam_id)

# SQANTI classification file- add sample names, calculate cpm on all transcripts before filtering
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

# Keep only annotated genes
sqanti_df_full = sqanti_df %>%
  filter(!is.na(associated_gene), !str_starts(associated_gene, "novelGene"))
#filter(!is.na(associated_gene), str_starts(associated_gene, "ENS"))

sqanti_df_full %<>% 
  left_join(gencode_gene, by = "associated_gene") %>% 
  left_join(gencode_df)

# combine with hashids, and selct cpm columns
hashids = read_tsv(mapping_file)

sqanti_ids = sqanti_df_full %>% 
  select(isoform_id = isoform, 
         ensg_gene_id = associated_gene, 
         any_of("gene_name"), 
         enst_transcript_id = associated_transcript, 
         any_of("transcript_name"), 
         any_of("gene_type"), 
         any_of("protein_id"), 
         ends_with("_counts"), 
         ends_with("_cpm"))

all_ids = hashids %>% select(isoform_id = transcript_id, hash_id) %>%
  inner_join(sqanti_ids, by = "isoform_id") 

write_tsv(all_ids, file.path(dir, paste0(basename, "_all_hashids_with_cpm.txt")))

# =============================================================================
# Apply filters sequentially
# =============================================================================

message("\n Starting filtering of SQANTI output...")

# Initialize dropout tracking
original_ids    = sqanti_df_full$isoform
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
  sqanti_df_full$dropout_reason[sqanti_df_full$isoform %in% dropped_ids] = "internal_priming"
  #dropout_tracker[["internal_priming"]] = dropped_ids
  
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
  sqanti_df_full$dropout_reason[sqanti_df_full$isoform %in% dropped_ids] = "template_switching"
  #dropout_tracker[["Template_switching_artifact"]] = dropped_ids
  
  message("Template switching filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts")
}

# Apply structural category filter
if (toupper(tclass_to_keep) == "ALL") {
  # Keep everything - no filtering
  ids      = sqanti_df$isoform
  kept_ids = sqanti_df$isoform
  message("No structural category filtering applied")
  
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
  
  message("Filtering to categories: ", paste(allowed_categories, collapse = ", "))
  
  sqanti_df %<>% 
    filter(structural_category %in% allowed_categories)
  
  kept_ids = sqanti_df$isoform
}
  
#allowed_categories = STRUCTURAL_CATEGORIES[[structural_level]]
dropped_ids = setdiff(ids, kept_ids)
sqanti_df_full$dropout_reason[sqanti_df_full$isoform %in% dropped_ids] = "structural_category"
#dropout_tracker[["structural_category_filtered"]] = dropped_ids
message("Structural category filter: kept ", length(kept_ids), " transcripts, dropped ", length(dropped_ids), " transcripts")

# Calculate final results
all_dropout_ids = sqanti_df_full$isoform[sqanti_df_full$dropout_reason != "kept"]

#all_dropout_ids = setdiff(original_ids, kept_ids)

#message(paste0("Final: ", length(kept_ids), " kept, ", length(all_dropout_ids), " dropped"))

# =============================================================================
# Write output files and generate summary
# =============================================================================

message("\n Writing output files...")

write_tsv(sqanti_df_full, file.path(dir, paste0(basename, "_classification_all.txt")))
write_tsv(sqanti_df, file.path(dir, paste0(basename, "_classification_filtered.txt")))

save_filtered_fasta(
  sqanti_fasta, 
  kept_ids, 
  file.path(dir, paste0(basename, "_corrected_filtered.fasta"))
)

# gtf and bed
# sample gtf
gtf = import(sqanti_gtf) %>%
  as.data.frame()

filtered_gtf = gtf[gtf$transcript_id %in% kept_ids, ]
gtf_ids = filtered_gtf %>% distinct(gene_id, transcript_id)

new_attributes = all_ids %>%
  filter(isoform_id %in% kept_ids) %>%
  select(transcript_id = isoform_id, everything()) %>%
  left_join(gtf_ids) %>%
  mutate(avg_cpm = rowMeans(across(contains("cpm")), na.rm = TRUE)) %>%
  group_by(gene_id) %>%
  mutate(
    gene_total_cpm = sum(avg_cpm, na.rm = TRUE),
    avg_ratio = round(avg_cpm / gene_total_cpm, 3)
  ) %>%
  ungroup() %>%
  select(transcript_id, reference_gene_name = gene_name, avg_ratio)

filtered_gtf %<>% 
  left_join(new_attributes, by = c("transcript_id")) %>%
  mutate(name = paste0(transcript_id, "|",
                       reference_gene_name, "|", 
                       avg_ratio))

gr_updated = makeGRangesFromDataFrame(filtered_gtf, keep.extra.columns = TRUE)
export(gr_updated, file.path(dir, paste0(basename, "_corrected_filtered.gtf")), format = "gtf")

# convert to bed12
gtf_to_bed12(gtf_path = file.path(dir, paste0(basename, "_corrected_filtered.gtf")),
             output_bed = file.path(dir, paste0(basename, "_corrected_filtered.bed")),
             color_by = "avg_ratio")

# Filtered hashid and cpm table
all_ids %>% filter(isoform_id %in% sqanti_df$isoform) %>%
  write_tsv(file.path(dir, paste0(basename, "_hashids_with_cpm_filtered.txt")))

# Print summary
cat("\n=== FILTERING SUMMARY ===\n")
cat(paste0("Sample: ", basename, "\n"))
cat(paste0("Original transcripts: ", length(original_ids), "\n"))
cat(paste0("Final transcripts: ", length(kept_ids), "\n"))
cat(paste0("Retention rate: ", round(100 * length(kept_ids) / length(original_ids), 1), "%\n"))