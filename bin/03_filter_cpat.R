#!/usr/bin/env Rscript

#' CPAT ORF Analysis Script - R Version
#' 
#' - Restores correct ID casing from CPAT's uppercased output using mapping file
#' - Maps ORF start/end to genomic coordinates via exon structure
#' - Checks for reference start codon matches
#' - Scores and classifies ORFs (Clear Best, Plausible, Low Quality, No Stop Codon)
#' - Generates GTF with CDS features for best ORFs 
#' 
#' Inputs:
#' - CPAT ORF fasta
#' - CPAT ORF probability tsv
#' - Filtered DNA fasta (from transcriptome subworkflow)
#' - Filtered GTF (from transcriptome subworkflow)
#' - Hashids mapping file (for restoring correct ID casing)
#'
#' Outputs:
#' - *.predicted_proteome.all_orfs_mapped.tsv: all ORFs with genomic coordinates and quality scores
#' - *.predicted_proteome.corrected_filtered_CDS.gtf: GTF with exon and CDS features for best ORF per transcript
#' 

# =============================================================================
# Load required libraries
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(Biostrings)
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
  make_option(c("--cpat_fasta"), type = "character", default = NULL,
              help = "Path to CPAT ORF sequences FASTA"),
  make_option(c("--cpat_results"), type = "character", default = NULL,
              help = "Path to CPAT ORF probability TSV"),
  make_option(c("--sample_fasta"), type = "character", default = NULL,
              help = "Path to corrected filtered DNA FASTA (transcriptome)"),
  make_option(c("--sample_gtf"), type = "character", default = NULL,
              help = "Path to corrected filtered GTF (transcriptome)"),
  make_option(c("--reference_gtf"), type = "character", default = NULL,
              help = "Path to reference GTF (GENCODE, RefSeq, etc.)"),
  make_option(c("--mapping_file"), type = "character", default = NULL,
              help = "Path to hashids mapping file for restoring correct ID casing"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Output directory for results"),
  make_option(c("--coding_threshold"), type = "double", default = 0.364,
              help = "CPAT coding probability threshold [default: %default]")
)

opt = parse_args(OptionParser(option_list = option_list))

required_args = c("basename", "cpat_fasta", "cpat_results", "sample_fasta",
                  "sample_gtf", "reference_gtf", "mapping_file", "output_dir")
missing = required_args[sapply(required_args, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste0("--", missing, collapse = ", "))
}

basename               = opt$basename
cpat_fasta_file        = opt$cpat_fasta
cpat_results_file      = opt$cpat_results
sample_full_fasta_path = opt$sample_fasta
sample_gtf_path        = opt$sample_gtf
gencode_gtf_path       = opt$reference_gtf
mapping_file           = opt$mapping_file
cpat_dir               = opt$output_dir
coding_threshold       = opt$coding_threshold

stopifnot("CPAT ORF FASTA file not found" = file.exists(cpat_fasta_file))
stopifnot("CPAT results file not found"   = file.exists(cpat_results_file))
stopifnot("Sample fasta file not found"   = file.exists(sample_full_fasta_path))
stopifnot("Sample gtf file not found"     = file.exists(sample_gtf_path))
stopifnot("Reference GTF not found"       = file.exists(gencode_gtf_path))
stopifnot("Mapping file not found"        = file.exists(mapping_file))

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
    ID       = tolower(str_split_fixed(names(orfs), "\t", 2)[,1]),
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
map_orfs_to_genome <- function(orf_coords, sample_exons, gencode_gtf, full_fasta) {
  
  # Count the number of upstream ORFs
  orf_coords %<>% left_join(full_fasta, by = "isoform_id") %>%
    mutate(upstream_atgs = ifelse(!is.na(full_sequence), 
                                  str_count(substr(full_sequence, 1, ORF_start - 1), "ATG"),
                                  Inf)
    ) %>%
    select(-full_sequence)
  
  # Map ORF Start coordinates to genomic positions- strandedness here is tricky
  mapped_starts = orf_coords %>%
    left_join(sample_exons, by = c("isoform_id"), relationship = "many-to-many") %>%
    filter(ORF_start > prior_length & ORF_start <= cumulative_length) %>% # which exon has the ORF start?
    mutate(orf_genomic_start = ifelse(strand == "+",
                                      start + (ORF_start - prior_length) - 1,
                                      end - (ORF_start - prior_length) + 1),
           orf_start_exon   = exon_number,
           orf_start_offset = ORF_start - prior_length  # position within this exon
    ) %>%
    select(-start, -end, -exon_length, -exon_number, -cumulative_length, -prior_length)
  
  # Map ORF End coordinates to genomic positions
  mapped_complete = mapped_starts %>%
    left_join(sample_exons, by = c("isoform_id", "seqnames", "strand"), relationship = "many-to-many") %>%
    filter(ORF_end > prior_length & ORF_end <= cumulative_length) %>% # which exon has the ORF end?
    mutate(
      orf_genomic_end = if_else(strand == "+",
                                start + (ORF_end - prior_length) - 1,
                                end - (ORF_end - prior_length) + 1),
      orf_end_exon   = exon_number,
      orf_end_offset = ORF_end - prior_length  # position within this exon
    ) %>%
    select(-start, -end, -exon_length, -exon_number, -cumulative_length, -prior_length)
  
  # Check for GENCODE start codon matches- looking for any match here
  gencode_start_codons = gencode_gtf %>%
    filter(type == "start_codon") %>% 
    mutate(gencode_pos = ifelse(strand == "+", start, end)) %>%
    select(seqnames, strand, gencode_pos) %>%
    distinct() %>%
    mutate(gencode_match = TRUE)
  
  mapped_orfs = mapped_complete %>% 
    left_join(gencode_start_codons, by = c("seqnames", "strand", "orf_genomic_start" = "gencode_pos")) %>%
    mutate(gencode_match = ifelse(is.na(gencode_match), FALSE, gencode_match))
  
  return(mapped_orfs)
}

#' Define plausible ORFs, can be more than one per transcript
#' From Plausible ORFs, select single best ORF per transcript
#' Logic: Prefer GENCODE matches with fewest upstream ATGs, otherwise highest coding score
#' @param mapped_orfs Data frame of mapped ORFs
#' @return Data frame with best ORF per transcript
call_best_orfs <- function(mapped_orfs) {
  
  # ORF quality scoring parameters but these don't seem to be used
  #coding_threshold  # CPAT protein-coding threshold
  #atg_shift  <- 10          # Sigmoid parameters for ATG penalty
  #atg_growth <- 0.5
  
  # per transcript, if there is a gencode match, pick the one with the fewest upstream ATGs
  mapped_orfs %<>% 
    mutate(#atg_penalty     = 1 - 1/(1 + exp(-atg_growth * (upstream_atgs - atg_shift))), 
           #composite_score = Coding_prob * (1 - atg_penalty),
           orf_quality     = case_when(
             Coding_prob <= coding_threshold ~ "Low Quality ORF",
             !has_stop_codon & !gencode_match ~ "No Stop Codon",
             TRUE ~ "Plausible ORF"
           ))
  
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

#' Get CDS coordinates by trimming exons to ORF boundaries
#' Uses the exon mapping information to extract only CDS regions
#' @param best_orfs Each row has an ORF with mapping info
#' @param sample_gtf Info on all exons
#' @return Data frame with CDS start/end coordinates
get_cds_coords <- function(best_orfs, sample_exons) {
  
  # Join mapped ORFs with their exons- select all CDS exons
  cds_all = best_orfs %>%
    left_join(sample_exons, by = c("isoform_id", "seqnames", "strand")) %>%
    filter(exon_number >= orf_start_exon & exon_number <= orf_end_exon) %>% # get all CDS exons
    group_by(isoform_id) %>%
    arrange(isoform_id, exon_number) %>%
    mutate(
      is_first_orf_exon = (exon_number == min(exon_number)),
      is_last_orf_exon = (exon_number == max(exon_number))
    ) %>%
    ungroup()
  
  # Trim first and last cds exons based on precise cds start and end
  cds_trimmed = cds_all %>%
    mutate(
      new_start = case_when(
        strand == "+" & is_first_orf_exon ~ start + orf_start_offset - 1,
        strand == "-" & is_last_orf_exon ~ end - orf_end_offset + 1,
        TRUE ~ start
      ),
      new_end = case_when(
        strand == "+" & is_last_orf_exon ~ start + orf_end_offset - 1,
        strand == "-" & is_first_orf_exon ~ end - orf_start_offset + 1,
        TRUE ~ end
      )
    ) %>%
    select(ID, isoform_id, seqnames, exon_number, new_start, new_end, strand)
  
  return(cds_trimmed)
}

#' Create updated sample GTF file that includes CDS type
#' Uses the exon mapping information to extract only CDS regions
#' @param sample_gtf Starting gtf file
#' @param all_cds_exons Trimmed CDS exon coordinates
#' @return GTF file of all transcripts with CDS features
write_gtf_with_cds <- function(sample_gtf, all_cds_exons) {
  
  # selecting type individually in case custom gtfs have more type columns- want to reset attributes column
  #sample_gtf %<>% 
    #select(isoform_id = transcript_id, everything())
    #left_join(gene_mapping %>% select(isoform_id, gencode_gene_id = gene_id, any_of("gene_name"), orf_isoform_id), by = "isoform_id")
  
  # Create transcript lines (one per transcript)
  transcript_lines = sample_gtf %>%
    filter(type == "transcript") %>%
    select(seqnames, type, start, end, strand, source, transcript_id, gene_id)
    #select(seqnames, type, start, end, strand, source, isoform_id, gencode_gene_id, any_of("gene_name"), orf_isoform_id)
  
  # Create exon lines
  exon_lines = sample_gtf %>%
    filter(type == "exon") %>%
    select(seqnames, type, start, end, strand, source, transcript_id, gene_id)
  
  # Create CDS lines
  cds_lines = all_cds_exons %>%
    #left_join(gene_mapping %>% select(isoform_id, gencode_gene_id = gene_id, any_of("gene_name"), orf_isoform_id), by = "isoform_id") %>%
    mutate(type = "CDS", source = sample_gtf$source[1]) %>%
    select(seqnames, type, start = new_start, end = new_end, strand, source, transcript_id = isoform_id, gene_id)
  
  # combine types and add attribute column
  updated_gtf = bind_rows(transcript_lines, exon_lines, cds_lines) %>%
    arrange(transcript_id, start) %>%
    mutate(
      attributes = if ("gene_name" %in% colnames(.)) {
        #paste0('gene_id "', gene_id, '"; transcript_id "', transcript_id, '"; gene_name "', gene_name, '"; ORF_id "', orf_isoform_id, '";')
        paste0('gene_id "', gene_id, '"; transcript_id "', transcript_id, '"; gene_name "', gene_name, '";')
      } else {
        paste0('gene_id "', gene_id, '"; transcript_id "', transcript_id, '";')
      },
      score = ".", 
      frame = "."
    ) %>%
    select(seqnames, source, type, start, end, score, strand, frame, attributes)
  
  cat(paste0("Combined ", nrow(transcript_lines), " transcript lines, ",
                 nrow(exon_lines), " exon lines, and ", 
                 nrow(cds_lines), " CDS lines"))
  
  return(updated_gtf)
}

# =============================================================================
# Main code broken down into steps
# =============================================================================

# === STEP 1: Load all input data ===
cat("\nSTEP 1: Loading input data")

# CPAT results, read mapping file for correct ID casing
id_lookup = read_tsv(mapping_file, show_col_types = FALSE) %>%
  select(isoform_id) %>%
  mutate(isoform_id_lower = tolower(isoform_id))

# CPAT results — tolower to fix CPAT's uppercasing, then restore correct casing
cpat_results = read_tsv(cpat_results_file, show_col_types = FALSE) %>%
  mutate(ID = tolower(ID)) %>%
  separate(ID, into = c("isoform_id_lower", "orf_rank"), sep = "_orf_", remove = FALSE) %>%
  left_join(id_lookup, by = "isoform_id_lower") %>%
  select(-isoform_id_lower)

# Load GTF files
sample_gtf  = import(sample_gtf_path) %>% as.data.frame()
gencode_gtf = import(gencode_gtf_path) %>% as.data.frame()

# Load sample fasta with full sequence- include UTRs, filtered transcripts from SQANTI
full = readDNAStringSet(sample_full_fasta_path)
full_fasta = tibble(
  isoform_id    = names(full),
  full_sequence = as.character(full)
)

# === STEP 2: Get dataframe of exons from gtf ===
cat("\nSTEP 2: Extracting exons from sample gtf")
sample_exons = sample_gtf %>%
  filter(type == "exon") %>%
  select(isoform_id = transcript_id, seqnames, start, end, strand) %>%
  mutate(exon_length = end - start + 1) %>%
  arrange(isoform_id, if_else(strand == "+", start, -start))

sample_exons %<>%
  group_by(isoform_id) %>%
  mutate(
    exon_number = row_number(),
    cumulative_length = cumsum(exon_length), # Calculate cumulative exon positions for coordinate mapping
    prior_length = lag(cumulative_length, default = 0)
  ) %>%
  ungroup()

# === STEP 3: Map ORFs to genomic coordinates ===
cat("\nSTEP 3: Mapping ORFs to genome")

stop_codon_status = check_stop_codons(cpat_fasta_file) # uses ORF (dna) fasta to return new column about stop codon 
orf_coords        = left_join(cpat_results, stop_codon_status, by = "ID")
mapped_orfs       = map_orfs_to_genome(orf_coords, sample_exons, gencode_gtf, full_fasta)

# === STEP 4: Define the best ORF per transcript ===
cat("\nSTEP 4: Calling best ORFs and filtering\n")

mapped_orfs_classified = call_best_orfs(mapped_orfs)
#mapped_orfs_classified %<>% left_join(gene_mapping, by = "isoform_id")

all_orfs       = mapped_orfs_classified
plausible_orfs = mapped_orfs_classified %>% filter(orf_quality == "Clear Best ORF" | orf_quality == "Plausible ORF")
best_orfs      = mapped_orfs_classified %>% filter(orf_quality == "Clear Best ORF")

#write_tsv(best_orfs, file.path(cpat_dir, paste0(basename, "_best_orfs_mapped.tsv")))
write_tsv(all_orfs, file.path(cpat_dir, paste0(basename, ".predicted_proteome.all_orfs_mapped.tsv")))

# === STEP 5: Write gtf for best ORFs, including CDS and exon types, no collapsing here ===
cat("\nSTEP 5: Writing GTF of best ORFs with CDS and exon types")
all_cds_exons = get_cds_coords(best_orfs, sample_exons)
# gene_mapping %<>% 
#   left_join(select(best_orfs, isoform_id, orf_isoform_id)) %>%
#   mutate(orf_isoform_id = ifelse(is.na(orf_isoform_id), "noORF", orf_isoform_id))

all_cds_exons %<>% 
  left_join(distinct(sample_gtf, transcript_id, gene_id), by = c("isoform_id" = "transcript_id"))
updated_gtf = write_gtf_with_cds(sample_gtf, all_cds_exons)

write.table(updated_gtf, file.path(cpat_dir, paste0(basename, ".predicted_proteome.corrected_filtered_CDS.gtf")), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

# =============================================================================
# Summary
# =============================================================================
all_sample_transcripts = distinct(sample_gtf, transcript_id, gene_id)

starting_transcripts = nrow(all_sample_transcripts)
n_transcripts_orfs   = n_distinct(all_orfs$isoform_id)
percent_orfs         = round((n_transcripts_orfs/starting_transcripts) * 100, 4)

n_cpat_orfs          = nrow(all_orfs)
n_plausible_orfs     = nrow(plausible_orfs)
percent_plausible    = round((n_plausible_orfs/n_cpat_orfs) * 100, 1)

n_best_orfs          = nrow(best_orfs %>% distinct(isoform_id))
n_noORFs             = starting_transcripts - n_best_orfs
percent_best         = round((n_best_orfs/starting_transcripts) * 100, 1)
percent_noORFs       = round((n_noORFs/starting_transcripts) * 100, 1)


cat("\n=== CPAT FILTER ANALYSIS COMPLETE ===\n")
cat(paste0("- CPAT called at least one ORF for ", n_transcripts_orfs, " / ", starting_transcripts, " (", percent_orfs,"%)", " transcripts.\n"))
cat(paste0("- Based on coding potential, ", n_plausible_orfs, " / ", n_cpat_orfs, " (", percent_plausible,"%)"," of CPAT ORFs were classified as plausible ORFs.\n")) 
cat(paste0("- ", n_noORFs, " / ", starting_transcripts, " (", percent_noORFs,"%)", " transcripts do not have a plausible ORF.\n")) 
cat(paste0("- After further filtering, a 'Clear Best ORF' was identified for ", n_best_orfs, " / ", starting_transcripts, " (", percent_best,"%)", " transcripts.\n")) 

cat("\n=== CPAT FILTER OUTPUT FILES ===\n")
cat(paste0("All CPAT ORFs with Quality Metrics: ", basename, ".predicted_proteome.all_cpat_orfs_mapped.tsv\n"))
#cat(paste0("The single best plausible ORF per transcript: ", basename, "_best_cpat_orfs_mapped.tsv"))
cat(paste0("GTF contains exon type for all transcripts and CDS type for transcripts with a best plausible ORF: ", basename, ".predicted_proteome.corrected_filtered_CDS.gtf"))
