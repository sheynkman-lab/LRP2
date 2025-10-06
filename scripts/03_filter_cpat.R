#!/usr/bin/env Rscript

#' ORF Analysis Script - R Version
#' 
#' Map ORF coordinates to genomic positions
#' Compares to GENCODE start codon annotations
#' Selects the best ORF for each transcript
#'
#' Outputs 

# =============================================================================
# Load required libraries
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(Biostrings)
  library(magrittr)
})

# =============================================================================
# Get environment variables
# =============================================================================

basename         <- Sys.getenv("OUTPUT_BASE_NAME")
sqanti_dir       <- file.path(Sys.getenv("OUTPUT_DIR"), "sqanti") 
cpat_dir         <- file.path(Sys.getenv("OUTPUT_DIR"), "orf_calling") 
gencode_gtf_path <- Sys.getenv("GENCODE_GTF_FILE")
coding_threshold <- as.numeric(Sys.getenv("CPAT_CODING_THRESHOLD", "0.364"))

# =============================================================================
# Check required files
# =============================================================================

cpat_fasta_file <- file.path(cpat_dir, paste0(basename, "_cpat.ORF_seqs.fa"))
stopifnot("Input ORF FASTA file not found" = file.exists(cpat_fasta_file))

cpat_results_file <- file.path(cpat_dir, paste0(basename, "_cpat.ORF_prob.tsv"))
stopifnot("CPAT results file not found" = file.exists(cpat_results_file))

sample_full_fasta_path <- file.path(sqanti_dir, paste0(basename, "_corrected_filtered.fasta"))
stopifnot("Sample full fasta file not found" = file.exists(sample_full_fasta_path))

sample_gtf_path <- file.path(sqanti_dir, paste0(basename, "_corrected_filtered.gtf"))
stopifnot("Sample gtf file not found" = file.exists(sample_gtf_path))

classification_file <- paste0(sqanti_dir, "/", basename, "_classification_filtered.txt")
stopifnot("File not found" = file.exists(classification_file))

# =============================================================================
# Helper functions
# =============================================================================

#' Check if ORF ends with stop codon
#' @param orf_fasta_file Path to ORF FASTA file
#' @return Data frame with ID and has_stop_codon columns
check_stop_codons <- function(orf_fasta_file) {
  
  stop_codons <- c("TAG", "TAA", "TGA")
  orfs        <- readDNAStringSet(orf_fasta_file)
  
  tibble(
    ID       = str_split_fixed(names(orfs), "\t", 2)[,1],
    sequence = as.character(orfs)
  ) %>%
    mutate(seq_length     = nchar(sequence),
           last_codon     = if_else(seq_length >= 3, str_sub(sequence, -3, -1), ""),
           has_stop_codon = last_codon %in% stop_codons) %>%
    select(ID, has_stop_codon)
}

#' Map ORF coordinates to genomic positions
#' This is a simplified version that focuses on the core mapping logic
#' @param orf_coords ORF coordinate data
#' @param sample_gtf GTF data for sample transcripts
#' @param gencode_gtf GTF data for GENCODE reference
#' @param sample_fasta Sample transcript sequences
#' @return Mapped ORF data with genomic coordinates
map_orfs_to_genome <- function(orf_coords, sample_gtf, gencode_gtf, full_fasta) {
  
  # Extract exon information from sample GTF
  sample_exons = sample_gtf %>%
    filter(type == "exon") %>%
    select(transcript_id, seqnames, start, end, strand) %>%
    mutate(exon_length = end - start + 1) %>%
    arrange(transcript_id, if_else(strand == "+", start, -start))
  
  # Calculate cumulative exon positions for coordinate mapping
  sample_exons %<>%
    group_by(transcript_id) %>%
    mutate(cumulative_length = cumsum(exon_length), prior_length = lag(cumulative_length, default = 0)) %>%
    ungroup()
  
  # Map ORF coordinates to genomic positions and pull number of upstream ATGs- strandedness here is tricky
  mapped_orfs = orf_coords %>%
    left_join(sample_exons, by = c("isoform_id" = "transcript_id")) %>%
    filter(ORF_start > prior_length & ORF_start <= cumulative_length) %>% # which exon has the ORF start?
    left_join(full_fasta, by = "isoform_id") %>%
    mutate(genomic_start = ifelse(strand == "+",
                                  start + (ORF_start - prior_length) - 1,
                                  end - (ORF_start - prior_length) + 1)) %>%
    mutate(upstream_atgs = ifelse(!is.na(full_sequence), 
                                  str_count(substr(full_sequence, 1, ORF_start - 1), "ATG"),
                                  Inf)) %>%
    select(-full_sequence)
  
  # Check for GENCODE start codon matches- looking for any match here
  gencode_start_codons = gencode_gtf %>%
    filter(type == "start_codon") %>% 
    mutate(gencode_pos = ifelse(strand == "+", start, end)) %>%
    select(seqnames, strand, gencode_pos) %>%
    distinct() %>%
    mutate(gencode_match = TRUE)
  
  mapped_orfs %<>% left_join(gencode_start_codons, by = c("seqnames", "strand", "genomic_start" = "gencode_pos")) %>%
    mutate(gencode_match = ifelse(is.na(gencode_match), FALSE, gencode_match))
  
  return(mapped_orfs)
}

#' Select best ORF for each transcript
#' Logic: Prefer GENCODE matches with fewest upstream ATGs, otherwise highest coding score
#' @param mapped_orfs Data frame of mapped ORFs
#' @return Data frame with best ORF per transcript
call_best_orfs <- function(mapped_orfs) {
  
  # ORF quality scoring parameters but these don't seem to be used
  coding_threshold  # CPAT protein-coding threshold
  atg_shift  <- 10          # Sigmoid parameters for ATG penalty
  atg_growth <- 0.5
  
  # per transcript, if there is a gencode match, pick the one with the fewest upstream ATGs
  mapped_orfs %<>% 
    mutate(atg_penalty     = 1 - 1/(1 + exp(-atg_growth * (upstream_atgs - atg_shift))), 
           composite_score = Coding_prob * (1 - atg_penalty),
           orf_quality     = ifelse(!has_stop_codon, "No Stop Codon", 
                                    ifelse(Coding_prob <= coding_threshold, "Low Quality ORF", "Plausible ORF")))
  
  best_orfs = mapped_orfs %>% 
    filter(orf_quality == "Plausible ORF") %>%
    group_by(isoform_id) %>% 
    arrange(desc(gencode_match), upstream_atgs, desc(Coding_prob)) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  mapped_orfs %<>%
    mutate(orf_quality = ifelse(ID %in% best_orfs$ID, "Clear Best ORF", orf_quality))
  
  return(mapped_orfs)
}

#' Translate ORF sequences to proteins and group by identical sequences
#' @param orfs Data frame with ORF coordinates
#' @param transcript_seqs Named vector of transcript sequences
#' @return Data frame with protein sequences and grouped transcript IDs
group_by_protein_sequence <- function(filtered_results, full_fasta) {
  
  orf_proteins = filtered_results %>%
    select(isoform_id, ORF_start, ORF_end) %>%
    left_join(full_fasta, by = "isoform_id") %>%
    filter(!is.na(full_sequence)) %>%
    mutate(orf_dna_sequence = substr(full_sequence, ORF_start, ORF_end))
  
  # Vectorized translation for speed
  dna      <- DNAStringSet(orf_proteins$orf_dna_sequence)
  proteins <- translate(dna, if.fuzzy.codon = "solve")
  
  orf_proteins %<>%
    mutate(orf_protein_sequence = as.character(proteins)) %>%
    filter(orf_protein_sequence != "") %>%
    select(isoform_id, orf_protein_sequence)
  
  # Group transcripts by identical protein sequences
  orf_groups = orf_proteins %>%
    group_by(orf_protein_sequence) %>%
    arrange(isoform_id) %>%
    mutate(
      orf_isoform_id = paste(isoform_id, collapse = " | "),
      orf_base_id    = isoform_id[1]  # First transcript as representative
    ) %>% 
    ungroup()
  
  return(orf_groups)
}

# =============================================================================
# Main code broken down into steps
# =============================================================================

# === STEP 1: Load all input data ===
message("\n--- STEP 1: Loading input data ---")

# CPAT results
cpat_results = read_tsv(cpat_results_file, show_col_types = FALSE) %>%
  separate(ID, into = c("isoform_id", "orf_rank"), sep = "_ORF_", remove = FALSE)

# Load GTF files
sample_gtf  = import(sample_gtf_path) %>% as.data.frame()
gencode_gtf = import(gencode_gtf_path) %>% as.data.frame()

# Load sample fasta with full sequence- include UTRs, filtered transcripts from SQANTI
full = readDNAStringSet(sample_full_fasta_path)
full_fasta = tibble(
  isoform_id    = names(full),
  full_sequence = as.character(full)
)

# SQANTI classification
classification = read_tsv(classification_file, show_col_types = FALSE)

# === STEP 2: Map ORFs to genomic coordinates ===
message("\n--- STEP 2: Mapping ORFs to genome ---")

stop_codon_status = check_stop_codons(cpat_fasta_file) # uses ORF (dna) fasta to return new column about stop codon 
orf_coords        = left_join(cpat_results, stop_codon_status, by = "ID")
mapped_orfs       = map_orfs_to_genome(orf_coords, sample_gtf, gencode_gtf, full_fasta)

# === STEP 3: Define the clear best ORF and print out all classified ORFs ===
message("\n--- STEP 3: Calling best ORFs and filtering---")

mapped_orfs_classified = call_best_orfs(mapped_orfs)
mapped_orfs_classified %<>% left_join(select(classification, 
                                             isoform_id = isoform, 
                                             structural_category, 
                                             gene_id = associated_gene, 
                                             transcript_id = associated_transcript,
                                             gene_name, 
                                             transcript_name))

original_n_transcripts = n_distinct(mapped_orfs_classified$isoform_id)

best_orfs = mapped_orfs_classified %>% 
  filter(orf_quality == "Clear Best ORF")

n_best_orfs = nrow(best_orfs %>% distinct(isoform_id))

write_tsv(mapped_orfs_classified, file.path(cpat_dir, paste0(basename, "_all_cpat_orfs_mapped.tsv")))

# === STEP 4: Group transcripts by ORF protein sequence ===
message("\n--- Grouping by ORF protein sequence ---")
orf_groups = group_by_protein_sequence(best_orfs, full_fasta)
best_orfs %<>% left_join(orf_groups)

write_tsv(best_orfs, file.path(cpat_dir, paste0(basename, "_best_cpat_orfs_mapped.tsv")))


# =============================================================================
# Summary
# =============================================================================

message("\n=== ANALYSIS COMPLETE ===")
message(paste0("Started with ORFs called for ", original_n_transcripts, " transcripts."))
message(paste0("After filtering, ", n_best_orfs, " transcripts have a 'Clear Best ORF'")) 
message(paste0("All mapped ORFs: ", basename, "_all_cpat_orfs_mapped.tsv"))
message(paste0("Best mapped ORFs: ", basename, "_best_cpat_orfs_mapped.tsv"))
message("Analysis completed successfully!")