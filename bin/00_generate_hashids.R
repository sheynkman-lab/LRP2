#!/usr/bin/env Rscript

#' Generate hashids for all transcripts from a gtf
#' 
#' - Convert gtf to psl format, 1-based to 0-based coordinate system
#' - Requires gene_id and transcript_id column in gtf/psl
#' - Hash id is calculated from psl coords (0-based)
#' - Mono-exonic transcripts do not have junction chain coords, so these are handled differently 
#' - Concatenate chr and strand to hashid
#'
#' Outputs mapping file with gene_id, transcript_id, and hash_id.

# =============================================================================
# Load required libraries
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(magrittr)
})

options(scipen = 999)

# =============================================================================
# Check required files
# =============================================================================

basename <- Sys.getenv("OUTPUT_BASE_NAME")
dir      <- file.path(Sys.getenv("OUTPUT_DIR"), "sqanti_transcript")
gtf      <- paste0(dir, "/", basename, "_corrected.gtf")

stopifnot("File not found" = file.exists(gtf))

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
  
  message("✓ Validation passed: No underscores in transcript_id or gene_id")
  
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

# =============================================================================
# Generating hashids
# =============================================================================

message("\n Generating hash ids for all transcripts...")

# Sample gtf to psl
psl = file.path(dir, paste0(basename, "_corrected.psl")) # output psl
convert_gtf_to_psl(gtf_input_path  = gtf, 
                   psl_output_file = psl)

# Run python script for hash ids- generates mapping file with transcript_id and hash_id
mapping_file = file.path(dir, paste0(basename, "_hashids_mapping.txt"))
hashlib_script = Sys.getenv("HASHLIB_SCRIPT")
system2("python", args = c(hashlib_script, psl, mapping_file))