#!/usr/bin/env Rscript

#' Generate hashids and Pipeline Transcript IDs
#' 
#' - Convert sample gtf to psl format (1-based to 0-based coordinate system)
#' - Requires gene_id and transcript_id column in gtf/psl (gene_name optional for readability)
#' - Hash id is calculated from junction coordinates via SHAKE-256 (0-based from psl, TSS/TES ignored)
#' - Mono-exonic transcripts do not have junction chain coords, handled differently (TSS/TES used)
#' - Concatenate chr and strand to hashid
#' - Assign new pipeline transcript IDs using SQANTI QC reference mapping and structural categories:
#'     FSM: GENE_NAME.ENST or GENE_NAME.NR (known isoform)
#'     Non-FSM: GENE_NAME.junction_hash (novel isoform)
#'     
#' Inputs:
#' - Sample gtf
#' - Reference gtf
#' - SQANTI classification file
#' 
#' Outputs:
#' - *_hashids_mapping.txt
#' 

# =============================================================================
# Load required libraries
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
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
  make_option(c("--sample_gtf"), type = "character", default = NULL,
              help = "Path to sample GTF (corrected, from SQANTI)"),
  make_option(c("--classification"), type = "character", default = NULL,
              help = "Path to SQANTI classification file"),
  make_option(c("--reference_gtf"), type = "character", default = NULL,
              help = "Path to reference GTF"),
  make_option(c("--hashlib_script"), type = "character", default = NULL,
              help = "Path to python hashlib_id_generator.py"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Output directory for results")
)

opt = parse_args(OptionParser(option_list = option_list))

required_args = c("basename", "sample_gtf", "classification", "reference_gtf", 
                  "hashlib_script", "output_dir")

missing = required_args[sapply(required_args, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste0("--", missing, collapse = ", "))
}

basename       = opt$basename
sample_gtf     = opt$sample_gtf
class_file     = opt$classification
reference_gtf    = opt$reference_gtf
hashlib_script = opt$hashlib_script
output_dir     = opt$output_dir

stopifnot("Sample GTF not found"          = file.exists(sample_gtf))
stopifnot("Classification file not found" = file.exists(class_file))
stopifnot("Reference GTF not found"         = file.exists(reference_gtf))
stopifnot("Hashlib script not found"      = file.exists(hashlib_script))

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
  
  # keep only the following columns 
  gtf_df %<>% 
    dplyr::select(seqnames, type, start, end, strand, transcript_id, gene_id)
  
  colnames(gtf_df) = c("chrom", "ty", "start", "end", "strand", "transcript_id", "gene_id")
  gtf_df$start     = gtf_df$start-1 # converts from 1-based to 0-based
  
  # Filter for exons only
  exons_df = gtf_df %>% 
    filter(ty == "exon")
  
  cols_to_convert           = c("chrom", "ty", "strand", "transcript_id", "gene_id")
  exons_df[cols_to_convert] = lapply(exons_df[cols_to_convert], as.character)
  
  # Check for underscores in gene_id and transcript_id
  if (any(grepl("_", exons_df$transcript_id, fixed = TRUE)) || 
      any(grepl("_", exons_df$gene_id, fixed = TRUE))) {
    stop("ERROR: Underscores '_' found in gene or transcript IDs. Remove from IDs in GTF.")
  }
  
  cat("\n✓ Validation passed: No underscores in transcript_id or gene_id")
  
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

#' Extract junction hash from full hash ID (remove chr and strand)
#' Full hash format: chr1_s3a2b1c4:e5f6g7h8_+
#' Returns: s3a2b1c4.e5f6g7h8
#' For monoexonic: chr1_monoexon:100-200_+ -> monoexon.100-200
extract_junction_hash = function(hash_id) {
  case_when(
    str_detect(hash_id, "monoexon") ~ str_replace(hash_id, "^[^_]+_(monoexon:[0-9]+-[0-9]+)_[+-]$", "\\1") %>%
      str_replace(":", "."),
    TRUE ~ str_replace(hash_id, "^[^_]+_(.+)_[+-]$", "\\1") %>%
      str_replace(":", ".")
  )
}

# =============================================================================
# Step 1: Generating hashids
# =============================================================================

cat("\nSTEP 1: Generating hash ids for all transcripts")

psl          = file.path(output_dir, paste0(basename, "_corrected.psl")) # output psl
hashid_file  = file.path(output_dir, paste0(basename, "_hashids_raw.txt"))

convert_gtf_to_psl(gtf_input_path  = sample_gtf, 
                   psl_output_file = psl)

# Run python script for hash ids- generates mapping file with transcript_id and hash_id
system2("python", args = c(hashlib_script, psl, hashid_file))

hashids = read_tsv(hashid_file, show_col_types = FALSE)
cat("\nGenerated hash IDs for ", nrow(hashids), " transcripts\n")

# =============================================================================
# Step 2: Read SQANTI classification and reference GTF
# =============================================================================

cat("\nSTEP 2: Reading SQANTI classification and GENCODE reference")

sqanti = read_tsv(class_file, show_col_types = FALSE) %>%
  filter(!is.na(associated_gene), !str_starts(associated_gene, "novelGene"))

# gene name lookup from reference
reference    = import(reference_gtf, format = "gtf")
reference_df = as.data.frame(reference)

reference_gene = reference_df %>%
  filter(type == "transcript") %>%
  select(associated_gene = gene_id,
         gene_name = any_of(c("gene_name", "gene")),
         gene_type = any_of(c("gene_type", "gene_biotype"))) %>%
  distinct()

reference_transcript = reference_df %>%
  filter(type == "transcript") %>%
  select(associated_gene = gene_id,
         associated_transcript = transcript_id,
         any_of("transcript_name")) %>%
  distinct()

# =============================================================================
# STEP 3: Build mapping table with new pipeline transcript IDs
# =============================================================================

cat("\nSTEP 3: Building pipeline transcript IDs")

# Combine hash IDs with SQANTI classification and gene names
mapping = hashids %>%
  dplyr::rename(original_transcript_id = transcript_id) %>%
  inner_join(sqanti %>% select(original_transcript_id = isoform,
                               associated_gene,
                               associated_transcript,
                               structural_category),
             by = "original_transcript_id") %>%
  left_join(reference_gene, by = "associated_gene") %>%
  left_join(reference_transcript, by = c("associated_gene", "associated_transcript"))

# Gene label: prefer gene_name, fall back to gene_id
mapping %<>%
  mutate(gene_label = if_else(!is.na(gene_name) & gene_name != "", gene_name, associated_gene))

# Junction hash (without chr and strand) for novel transcript IDs
mapping %<>%
  mutate(junction_hash = extract_junction_hash(hash_id))

# New transcript ID:
#   FSM -> GENE_NAME.reference_transcript_id
#   Non-FSM -> GENE_NAME.junction_hash
mapping %<>%
  mutate(
    isoform_id = case_when(
      structural_category == "full-splice_match" & !is.na(associated_transcript) ~
        paste0(gene_label, "::", associated_transcript),
      TRUE ~
        paste0(gene_label, "::", junction_hash)
    )
  )

# =============================================================================
# STEP 4: Check for redundant isoform IDs
# =============================================================================

cat("\nSTEP 4: Checking for isoform ID redundancy")

# Extract chromosome from hash_id for chrX/Y redundancy
mapping %<>% mutate(chrom = str_extract(hash_id, "^[^_]+"))

dupe_ids = mapping %>%
  group_by(isoform_id) %>%
  filter(n() > 1) %>%
  ungroup()

if (nrow(dupe_ids) > 0) {
  
  # same isoform_id but different chromosomes (e.g., chrX/chrY)
  chr_dupes = dupe_ids %>%
    group_by(isoform_id) %>%
    filter(n_distinct(chrom) > 1) %>%
    ungroup()
  
  # same isoform_id and same chromosome (different TSS/TES)
  other_dupes = dupe_ids %>%
    group_by(isoform_id) %>%
    filter(n_distinct(chrom) == 1) %>%
    ungroup()
  
  if (nrow(chr_dupes) > 0) {
    n_chr = n_distinct(chr_dupes$isoform_id)
    cat("Chr gene redundancy e.g., chrX/chrY: ", n_chr, " isoform IDs on multiple chromosomes (",
            nrow(chr_dupes), " total rows). Appending chromosome to gene label.")
    
    chr_fix = chr_dupes %>%
      mutate(isoform_id = paste0(gene_label, ".", chrom, "::", junction_hash))
    
    # chr_fix = chr_dupes %>%
    #   mutate(isoform_id = str_replace(isoform_id, "^([^:]+)", paste0("\\1(", chrom, ")")))
    
    mapping = mapping %>%
      filter(!original_transcript_id %in% chr_dupes$original_transcript_id) %>%
      bind_rows(chr_fix)
  }
  
  if (nrow(other_dupes) > 0) {
    n_other = n_distinct(other_dupes$isoform_id)
    cat("TSS/TES redundancy: ", n_other, " isoform IDs with same junctions in the same gene, different ends (",
            nrow(other_dupes), " total rows). Appending ::1, ::2, etc.")
    
    other_fix = other_dupes %>%
      group_by(isoform_id) %>%
      mutate(rn = row_number()) %>%
      ungroup() %>%
      mutate(isoform_id = paste0(isoform_id, "::", rn)) %>%
      select(-rn)
    
    mapping = mapping %>%
      filter(!original_transcript_id %in% other_dupes$original_transcript_id) %>%
      bind_rows(other_fix)
  }
  
} else {
  cat("No isoform ID redundancy detected")
}

# =============================================================================
# STEP 5: Write output mapping file
# =============================================================================

cat("\nSTEP 5: Writing mapping file\n")

output_mapping = mapping %>%
  mutate(
    reference_transcript_id = if_else(structural_category == "full-splice_match",
                                      associated_transcript, NA_character_)
  )

if ("transcript_name" %in% colnames(mapping)) {
  output_mapping %<>%
    mutate(transcript_name = if_else(structural_category == "full-splice_match",
                                     transcript_name, NA_character_))
}

output_mapping %<>%
  select(isoform_id,
         original_transcript_id,
         hash_id,
         reference_gene_id = associated_gene,
         any_of("gene_name"),
         any_of("gene_type"),
         reference_transcript_id,
         any_of("transcript_name"),
         structural_category)

mapping_output = file.path(output_dir, paste0(basename, "_hashids_mapping.txt"))
write_tsv(output_mapping, mapping_output)

n_total = nrow(output_mapping)
n_fsm   = sum(output_mapping$structural_category == "full-splice_match")
n_novel = n_total - n_fsm
n_genes = n_distinct(output_mapping$reference_gene_id)

cat("\n=== HASH ID GENERATION COMPLETE ===\n")
cat("Total transcripts: ", n_total, " across ", n_genes, " genes\n")
cat("FSM (reference-based IDs): ", n_fsm, "\n")
cat("Non-FSM (hash-based IDs): ", n_novel, "\n")
cat("Output: ", mapping_output, "\n")
