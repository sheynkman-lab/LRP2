#!/usr/bin/env Rscript

# This script inputs gencode gtf and fasta files to create a unified reference table

# Load required libraries
library(tidyverse)
library(rtracklayer)
library(Biostrings)
library(argparser)

# ==============================================================================
# COMMAND LINE ARGUMENT SETUP
# ==============================================================================

p <- arg_parser("Process GTF and FASTA files to create reference tables")

# Input files arguments
p <- add_argument(p, "--gtf", help = "Gencode GTF input file location")
p <- add_argument(p, "--fa", help = "Gencode FASTA input file location") 

args <- parse_args(p)

# Check if files exist
if (!file.exists(args$gtf)) stop("GTF file not found: ", args$gtf)
if (!file.exists(args$fa)) stop("FASTA file not found: ", args$fa)


# ==============================================================================
# FUNCTION 1: GENERATE GENE/TRANSCRIPT MAPPING TABLE FROM GTF
# ==============================================================================

# Import GTF file using rtracklayer
# This handles the complex GTF format and extracts all attributes properly

# temp
gtf_file = "/home/cwp5au/sheynkman/external_data/GENCODE_v48/gencode.v48.basic.annotation.gtf"
  
gene_transcript_map <- function(gtf_file, args) {
  
  # Import GTF file
  gtf <- import(gtf_file)
  gtf_df <- as.data.frame(gtf)
  
  # Filter for transcripts only
  transcripts <- gtf_df %>%
    filter(type == "transcript") %>%
    select(gene_id, gene_name, transcript_id, transcript_name, transcript_type, protein_id)
  
  # Write output files
  # write_tsv(transcripts, args$ensg_gene, col_names = FALSE)
  
}

# ==============================================================================
# FUNCTION 2: CREATE TRANSCRIPT AND GENE LENGTH TABLE FROM FASTA
# ==============================================================================

iso_len_tab <- function(fa_file, args) {
  message("Processing FASTA file...")
  
  # Read FASTA file
  fasta <- readDNAStringSet(fa_file)
  
  # Extract information from headers
  headers <- names(fasta)
  
  iso_data <- tibble(header = headers) %>%
    separate(header, into = paste0("col", 1:10), sep = "\\|", fill = "right") %>%
    transmute(
      isoform = str_remove(col5, '"'),
      gene = str_remove(col6, '"'),
      length = as.numeric(str_remove(col7, '"'))
    ) %>%
    filter(!is.na(isoform), !is.na(gene), !is.na(length))
  
  write_tsv(iso_data, args$isoname_lens)
  message("Isoform length table created")
  
  return(iso_data)
}

gene_len_tab <- function(iso_data, args) {
  message("Creating gene length statistics...")
  
  gene_stats <- iso_data %>%
    group_by(gene) %>%
    summarise(
      avg_len = round(mean(length), 1),
      min_len = min(length),
      max_len = max(length),
      .groups = "drop"
    )
  
  write_tsv(gene_stats, args$gene_lens)
  message("Gene length statistics table created")
}

# ==============================================================================
# FUNCTION 4: EXTRACT PROTEIN-CODING GENES LIST
# ==============================================================================

protein_coding_genes <- function(gtf_file, args) {
  message("Extracting protein coding genes...")
  
  gtf <- import(gtf_file)
  gtf_df <- as.data.frame(gtf)
  
  pc_genes <- gtf_df %>%
    filter(gene_type == "protein_coding") %>%
    pull(gene_name) %>%
    unique() %>%
    sort()
  
  writeLines(pc_genes, args$protein_coding_genes)
  message("Protein coding genes list created")
}

# ==============================================================================
# MAIN EXECUTION FUNCTION
# ==============================================================================

main <- function() {
  message("Starting reference table generation...")
  message("Input GTF file: ", args$gtf)
  message("Input FASTA file: ", args$fa)
  message("=" %>% str_dup(50))
  
  # Step 1: Create gene/transcript mapping tables from GTF
  gen_map(args$gtf, args)
  
  # Step 2: Process FASTA file to create transcript length table
  iso_data <- iso_len_tab(args$fa, args)
  
  # Step 3: Generate gene-level length statistics from transcript data
  gene_len_tab(iso_data, args)
  
  # Step 4: Extract list of protein-coding genes
  protein_coding_genes(args$gtf, args)
  
  message("=" %>% str_dup(50))
  message("All reference tables created successfully!")
  message("Output files:")
  message("  - ENSG->Gene: ", args$ensg_gene)
  message("  - ENST->Isoname: ", args$enst_isoname) 
  message("  - Gene->ENSP: ", args$gene_ensp)
  message("  - Gene->Isoname: ", args$gene_isoname)
  message("  - Transcript lengths: ", args$isoname_lens)
  message("  - Gene length stats: ", args$gene_lens)
  message("  - Protein coding genes: ", args$protein_coding_genes)
}