#!/usr/bin/env Rscript

#' Protein Classification and Filtering
#' 
#' Classifies and filters long-read derived protein isoforms by integrating
#' SQANTI-protein output with 5' UTR structure analysis and ORF predictions.
#' 
#' STEP 1: Process GENCODE reference
#' Build exon chain strings and GRanges for overlap/junction matching
#'   
#' STEP 2: Process long-read transcripts from CDS GTF
#' Extract 5' UTR structure, CDS boundaries, stop codon positions,
#' and distance from stop codon to first downstream junction
#'   
#' STEP 3: Classify 5' UTRs
#'  - "subset": 5' UTR is contained within annotated boundaries
#'    (monoexonic TSS within known exon, or multiexonic UTR junctions
#'     match reference and TSS does not extend beyond the matched 
#'     reference start position)
#'  - "unique": 5' UTR extends beyond annotated boundaries 
#'    (monoexonic TSS protrudes beyond known exon, or multiexonic UTR 
#'     has novel junction chain or upstream extension)
#'     
#' STEP 4: Protein classification
#' Integrates SQANTI-protein N-term/splice/C-term comparisons with 
#' 5' UTR classification to assign each isoform to:
#'  - FPM (Full Protein Match): exact match to reference protein
#'  - IPM (Incomplete Protein Match): N-terminal truncation
#'  - NPC (Novel Protein Combination): known elements in novel combination
#'  - NPE (Novel Protein Element): at least one novel N-term, splice, or C-term
#'   
#' STEP 5: Reconcile gene IDs and names
#' Resolve discrepancies between transcript-level and protein-level gene assignments
#' 
#' STEP 6: Filter to high-confidence protein isoforms
#' Filters by: no ORF, NMD (with rescue for single junction within 
#' threshold distance of stop codon), atypical SQANTI categories, 
#' and protein classification
#' 
#' STEP 7: Write ORF-centric outputs
#' Group transcripts by identical protein sequence, write high-confidence
#' ORF GTF, BED12 (with viridis coloring by expression ratio), 
#' CPM count table, and full proteome FASTA
#'
#' Inputs: 
#'  - GENCODE GTF
#'  - sample CDS GTF 
#'  - DNA FASTA
#'  - mapped ORFs
#'  - SQANTI protein classification
#'  - CPM counts
#'  
#' Outputs: 
#' One row per transcript (not collapsed)
#'  - *.predicted_proteome.best_ORF_summary.txt
#'  - *.predicted_proteome.best_ORF.fa (for proteome reference building)
#' 
#' One row per unique ORF (collapsed)
#'  - *.predicted_proteome.collapsed_high_confidence_ORF_hashids_with_cpm.txt
#'  - *.predicted_proteome.collapsed_high_confidence_ORF.gtf
#'  - *.predicted_proteome.collapsed_high_confidence_ORF.bed
#'


# =============================================================================
# Load libraries
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(Biostrings)
  library(tidyverse)
  library(magrittr)
  library(rtracklayer)
  library(viridis)
  library(data.table)
  library(optparse)
})

options(scipen = 999)

# =============================================================================
# Get environment variables and check required files
# =============================================================================

option_list = list(
  make_option(c("--basename"), type = "character", default = NULL,
              help = "Output base name"),
  make_option(c("--gencode_gtf"), type = "character", default = NULL,
              help = "Path to GENCODE GTF file"),
  make_option(c("--sample_cds_gtf"), type = "character", default = NULL,
              help = "Path to sample CDS GTF file from orf calling module"),
  make_option(c("--sample_dna_fasta"), type = "character", default = NULL,
              help = "Path to corrected filtered DNA FASTA from transcriptome subworkflow"),
  make_option(c("--mapped_orfs"), type = "character", default = NULL,
              help = "Path to all ORFs mapped TSV"),
  make_option(c("--protein_sqanti"), type = "character", default = NULL,
              help = "Path to SQANTI protein classification TSV"),
  make_option(c("--cpm_file"), type = "character", default = NULL,
              help = "Path to hashids with CPM filtered file from transcriptome subworkflow"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Output directory for results"),
  make_option(c("--min_junctions_after_stop"), type = "integer", default = 0,
              help = "Minimum junctions after stop codon for NMD filter [default: %default]"),
  make_option(c("--protein_class_keep"), type = "character", default = "FPM,NPC,NPE",
              help = "Comma-separated protein classes to keep [default: %default]"),
  make_option(c("--nmd_rescue_dist"), type = "integer", default = 25,
              help = "Max distance (bp) from stop to junction for NMD rescue [default: %default]")
)

opt = parse_args(OptionParser(option_list = option_list))

required_args = c("basename", "gencode_gtf", "sample_cds_gtf", "sample_dna_fasta", 
                  "mapped_orfs", "protein_sqanti", "cpm_file", "output_dir")

missing = required_args[sapply(required_args, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste0("--", missing, collapse = ", "))
}

basename                       = opt$basename
gencode_gtf_path               = opt$gencode_gtf
sample_cds_gtf_path            = opt$sample_cds_gtf
sample_dna_fasta_path          = opt$sample_dna_fasta
mapped_orfs                    = opt$mapped_orfs
protein_sqanti_path            = opt$protein_sqanti
cpm_file_path                  = opt$cpm_file
output_dir                     = opt$output_dir
min_junctions_after_stop_codon = opt$min_junctions_after_stop
nmd_rescue_dist                = opt$nmd_rescue_dist
pclass_base_to_keep            = strsplit(opt$protein_class_keep, ",")[[1]]


# =============================================================================
# Helper Functions
# =============================================================================

#' Classify multiexonic 5' UTR based on junction chain
#' @param tss Integer: TSS coordinate
#' @param strand Character: '+' or '-'
#' @param junc_chain Character: junction chain string
#' @param gc_chains Character vector: GENCODE junction chains for this gene
#' @return Character: "perfect_subset", "known_protruding", or "novel"
classify_multiexonic_utr = function(tss, strand, junc_chain, gc_chains) {
  
  # Check for NA or invalid inputs
  if (is.na(tss) || is.na(strand) || is.na(junc_chain) || is.na(gc_chains)) return("novel")
  if (junc_chain == "" || length(gc_chains) == 0) return("novel")
  
  # Split the concatenated chains
  gc_chains = str_split(gc_chains, "\\|\\|\\|")[[1]]
  
  status = "novel"
  
  for (gc_chain in gc_chains) {
    
    # Check for NA values in gc_chain
    if (is.na(gc_chain)) next
    
    # Check if junction chain is detected
    match_result = str_detect(gc_chain, fixed(junc_chain))
    if (is.na(match_result) || !match_result) next
    
    # Junction chain matches - check if TSS protrudes beyond reference
    if (strand == "+") {
      prefix  = str_split(gc_chain, fixed(junc_chain))[[1]][1]
      coords  = str_split(prefix, "_|-")[[1]][1]
      gc_5end = as.integer(coords[length(coords)])
      
      if (is.na(gc_5end)) next
      
      if ((tss + 9) < gc_5end) {
        status = "known_protruding"
      } else {
        return("perfect_subset")
      }
    } else {
      suffix  = str_split(gc_chain, fixed(junc_chain))[[1]][2]
      coords  = str_split(suffix, "_|-")[[1]][1]
      gc_5end = as.integer(coords[2])
      
      if (is.na(gc_5end)) next
      
      if ((tss - 9) > gc_5end) {
        status = "known_protruding"
      } else {
        return("perfect_subset")
      }
    }
  }
  
  return(status)
}

#' Translate ORF sequences to proteins and group by identical sequences
#' @param orf_info Data frame with ORF coordinates, classification, counts and cpm
#' @param full_fasta Named vector of transcript sequences
#' @return Data frame with protein sequences and grouped transcript IDs
group_by_protein_sequence <- function(orf_info, full_fasta) {
  
  orf_proteins = orf_info %>%
    left_join(full_fasta, by = "transcript_id") %>%
    filter(!is.na(full_dna_sequence)) %>%
    mutate(orf_dna_sequence = substr(full_dna_sequence, ORF_start, ORF_end))
  
  # Vectorized translation for speed
  dna      <- DNAStringSet(orf_proteins$orf_dna_sequence)
  proteins <- translate(dna, if.fuzzy.codon = "solve")
  
  orf_proteins %<>%
    mutate(orf_aa_sequence = as.character(proteins)) %>%
    filter(orf_aa_sequence != "") %>%
    select(-orf_dna_sequence, -full_dna_sequence)
  
  # Group transcripts by identical protein sequences, base id should be high confidence isoform with the highest expression
  orf_groups = orf_proteins %>%
    mutate(avg_cpm = rowMeans(select(., contains("cpm")), na.rm = TRUE)) %>%
    mutate(filter_status = factor(filter_status, levels = c("high_confidence", "NMD", "sqanti_classification", "sqanti_atypical"))) %>%
    group_by(orf_aa_sequence, gene_id) %>%
    arrange(filter_status, desc(avg_cpm), transcript_id, .by_group = TRUE) %>%
    mutate(
      orf_all_isoform_id = paste(transcript_id, collapse = ","),
      #orf_hc_isoform_id  = na_if(paste(transcript_id[filter_status == "high_confidence"], collapse = ","), ""),
      orf_base_id        = transcript_id[1]  # First transcript as representative
    ) %>% 
    ungroup()
  
  return(orf_groups)
}

#' Convert GTF to colored BED12 for genome browser visualization
#'
#' @param gtf_path Path to input GTF file
#' @param output_bed Path to output BED12 file
#' @param color_by Column name for coloring (numeric, 0-1). High values = purple, low = yellow. Default: NULL (black)
gtf_to_bed12 <- function(gtf_path, output_bed, color_by = NULL, track_name = NULL) {
  
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
  } else {
    ex <- dt[type == "CDS"]       # Use CDS for blocks (CDS only)
  }
  
  cds <- dt[type == "CDS"]        # Always get CDS for thick regions
  
  # Set keys for fast joins
  setkey(ex, orf_base_id)
  if(nrow(cds) > 0) setkey(cds, orf_base_id)
  
  # Sort transcripts by gene, then by ratio (descending)
  if(!is.null(color_by)) {
    tx[, sort_val := as.numeric(get(color_by))]
    tx <- tx[order(gene_id, -sort_val, orf_base_id)]
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
  
  tx$name <- gsub(" ", "_", tx$name)
  
  # Build BED12 rows
  bed_list <- vector("list", nrow(tx))
  
  for(i in 1:nrow(tx)) {
    t <- tx[i]
    e <- ex[.(t$orf_base_id)][order(start)]
    
    if(nrow(e) == 0) next
    
    # Get CDS boundaries for thick regions
    if(nrow(cds) > 0) {
      t_cds <- cds[.(t$orf_base_id)]
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
  
  # Write track header + BED data
  if (!is.null(track_name)) {
    track_line = paste0('track name="', track_name, '" ',
                        'description="', track_name, '" ',
                        'itemRgb=On')
    writeLines(track_line, output_bed)
  }
  
  # Remove NULLs and write
  bed <- do.call(rbind, bed_list[!sapply(bed_list, is.null)])
  write.table(bed, output_bed, sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
  
}

# =============================================================================
# STEP 1: Process GENCODE reference
# =============================================================================

cat("\n--- STEP 1: Processing GENCODE reference ---\n")

# Load GENCODE GTF
gencode_gtf = import(gencode_gtf_path) %>% as.data.frame()

# Identify protein-coding transcripts (have CDS)
pc_transcripts = gencode_gtf %>%
  filter(type == "CDS") %>%
  pull(transcript_id) %>%
  unique()

# Extract exons and create GRanges for overlap checking
gencode_exons = gencode_gtf %>%
  filter(type == "exon", transcript_id %in% pc_transcripts) %>%
  select(seqnames, start, end, gene_id, transcript_id, strand)

gencode_gr = gencode_exons %>%
  select(seqnames, start, end, strand) %>%
  makeGRangesFromDataFrame()

gencode_gr = GenomicRanges::reduce(gencode_gr) # Reduces overlapping exons into merged genomic intervals

# Create exon chain strings per gene for junction matching (could potentially be tweaked to mimic collapse logic)
gencode_chains = gencode_exons %>%
  arrange(gene_id, transcript_id, start) %>%
  group_by(gene_id, transcript_id) %>%
  summarise(
    exon_chain = paste(paste0(start, "-", end), collapse = "_"),
    .groups = "drop_last"  # Keep gene grouping!
  ) %>%
  summarise(chains = paste(exon_chain, collapse = "|||"), .groups = "drop") 

# =============================================================================
# STEP 2: Process long read transcripts
# =============================================================================

cat("\n--- STEP 2: Processing Sample transcripts ---\n")

sample_gtf = import(sample_cds_gtf_path) %>% as.data.frame()

# Create TxDb and extract 5' UTRs directly
txdb      = makeTxDbFromGFF(sample_cds_gtf_path, format = "gtf")
utr_by_tx = fiveUTRsByTranscript(txdb, use.names = TRUE)

# Get transcript metadata (gene name, chromosome, strand)
tx_metadata = sample_gtf %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id, seqnames, strand)

# Get CDS Nterm and TSS info 
sample_cds = sample_gtf %>%
  filter(type == "CDS") %>%
  group_by(transcript_id) %>%
  summarise(
    nterm = ifelse(first(strand) == "+", min(start), max(end)),
    stop_codon = ifelse(first(strand) == "+", max(end), min(start)),
    .groups = "drop"
  )

sample_exons = sample_gtf %>%
  filter(type == "exon") %>%
  group_by(transcript_id) %>%
  summarise(
    tss = ifelse(first(strand) == "+", min(start), max(end)),
    .groups = "drop"
  )

# Get junctions (intron boundaries) from exons and compute distance from stop codon to first downstream junction
exon_junctions = sample_gtf %>%
  filter(type == "exon") %>%
  group_by(transcript_id) %>%
  arrange(start, .by_group = TRUE) %>%
  summarise(
    strand = first(strand),
    # Junction positions are at exon-exon boundaries (end of exon_i, start of exon_i+1)
    junc_positions = list(data.frame(
      junc_left = end[-n()],
      junc_right = start[-1]
    )),
    .groups = "drop"
  ) %>%
  unnest(junc_positions) %>%
  left_join(sample_cds, by = "transcript_id") %>%
  mutate(
    # Distance from stop codon to the nearest junction downstream of stop
    dist_stop_to_junc = case_when(
      strand == "+" & junc_left >= stop_codon ~ junc_left - stop_codon,
      strand == "-" & junc_right <= stop_codon ~ stop_codon - junc_right,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(dist_stop_to_junc)) %>%
  group_by(transcript_id) %>%
  summarise(
    dist_stop_to_first_junc = min(dist_stop_to_junc),
    .groups = "drop"
  )

# Convert 5' UTR GRangesList to tibble with stats
utr_stats = tibble(
  transcript_id    = names(utr_by_tx),
  num_5utr_exons   = lengths(utr_by_tx),
  utr_length       = sum(width(utr_by_tx))
)

# Create chain strings for 5' UTRs
utr_chains = tibble(
  transcript_id = names(utr_by_tx),
  utr_chain = sapply(utr_by_tx, function(gr) {
    if (length(gr) == 0) return("")
    gr = sort(gr)
    paste(paste0(start(gr), "-", end(gr)), collapse = "_")
  })
) %>%
  mutate(
    junc_chain = case_when(
      utr_chain == "" ~ "",
      !str_detect(utr_chain, "_") ~ "",
      TRUE ~ {
        coords = str_extract_all(utr_chain, "\\d+")
        sapply(coords, function(x) {
          if (length(x) <= 2) return("")
          junctions = as.numeric(x[2:(length(x)-1)])
          paste(junctions, collapse = "_")
        })
      }
    )
  )

# Combine all info
utr_info = tx_metadata %>%
  left_join(sample_exons, by = "transcript_id") %>%
  left_join(sample_cds, by = "transcript_id") %>%
  left_join(utr_stats, by = "transcript_id") %>%
  left_join(utr_chains, by = "transcript_id") %>%
  mutate(
    # Handle transcripts with no UTR
    num_5utr_exons = replace_na(num_5utr_exons, 0),
    utr_length = replace_na(utr_length, 0),
    utr_exon_status = case_when(
      num_5utr_exons == 0 ~ "no_utr",
      num_5utr_exons == 1 ~ "monoexonic",
      TRUE ~ "multiexonic"
    )
  )

# =============================================================================
# STEP 3: Classify 5' UTRs and merge with SQANTI protein
# =============================================================================

cat("\n--- STEP 3: Classifying 5' UTRs ---\n")

# Read in SQANTI protein classification to get gene map to reference
sqanti_class = read_tsv(protein_sqanti_path, show_col_types = FALSE)
gene_map = sqanti_class %>% 
  select(transcript_id = isoform_id, reference_gene_id = tx_gene) %>%
  distinct()

# Check TSS overlap with GENCODE using GRanges 
tss_gr = GRanges(
  seqnames = utr_info$seqnames,
  ranges   = IRanges(start = utr_info$tss, end = utr_info$tss),
  strand   = utr_info$strand
)

overlaps = findOverlaps(tss_gr, gencode_gr)
utr_info$tss_in_gc_exons = seq_len(nrow(utr_info)) %in% queryHits(overlaps)

# For monoexonic: calculate distance to nearest GENCODE exon 
distances       = rep(NA_real_, nrow(utr_info))
non_overlapping = which(!utr_info$tss_in_gc_exons & utr_info$utr_exon_status == "monoexonic")

if (length(non_overlapping) > 0) {
  nearest                    = distanceToNearest(tss_gr[non_overlapping], gencode_gr)
  distances[non_overlapping] = mcols(nearest)$distance
}

utr_info$distance_to_gc = distances

# Join GENCODE chains for multiexonic comparison
utr_info = utr_info %>%
  left_join(gene_map, by = "transcript_id") %>%
  left_join(select(gencode_chains, reference_gene_id = gene_id, everything()), by = "reference_gene_id")

# Classify monoexonic UTRs
utr_info = utr_info %>%
  mutate(
    mono_class = case_when(
      utr_exon_status != "monoexonic" ~ NA_character_,
      tss_in_gc_exons ~ "within",
      distance_to_gc >= 10 ~ "protruding",
      TRUE ~ "within"
    )
  )

# Classify multiexonic UTRs (needs rowwise for junction matching)
multi_idx = which(utr_info$utr_exon_status == "multiexonic")

if (length(multi_idx) > 0) {
  multi_results = utr_info[multi_idx, ] %>%
    rowwise() %>%
    mutate(
      multi_class = classify_multiexonic_utr(tss, strand, junc_chain, chains)
    ) %>%
    ungroup() %>%
    pull(multi_class)
  
  utr_info$multi_class            = NA_character_
  utr_info$multi_class[multi_idx] = multi_results
}

# Final classification 
utr_results = utr_info %>%
  mutate(
    junc_cat = case_when(
      is.na(nterm) ~ "no_orf",
      utr_exon_status == "no_utr" ~ "no_utr",
      utr_exon_status == "monoexonic" ~ mono_class,
      utr_exon_status == "multiexonic" ~ multi_class,
      TRUE ~ "unknown"
    ),
    
    # Map to final categories
    utr_cat = case_when(
      junc_cat %in% c("no_orf", "no_utr", "unknown") ~ junc_cat,
      junc_cat %in% c("within", "perfect_subset") ~ "subset",
      junc_cat %in% c("protruding", "known_protruding", "novel") ~ "unique",
      TRUE ~ "unknown"
    )
  ) %>%
  select(transcript_id, num_5utr_exons, utr_exon_status, tss_in_gc_exons, junc_cat, utr_cat)

# Merge with SQANTI protein, adjust C term differences since CPAT includes stop codon in coords and gencode does not
sqanti_class %<>% mutate(
  pr_cterm_diff = ifelse(!is.na(pr_cterm_diff), pr_cterm_diff + 3, pr_cterm_diff),
  pr_cterm_gene_diff = ifelse(!is.na(pr_cterm_gene_diff), pr_cterm_gene_diff - 3, pr_cterm_gene_diff),
  pr_chang = ifelse(!is.na(pr_chang), pr_chang - 3, pr_chang))

utr_output = sqanti_class %>%
  select(transcript_id = isoform_id, everything()) %>%
  full_join(utr_results, by = "transcript_id")

# 5' UTR summary statistics
cat("\n--- 5'UTR Classification complete! ---\n")
cat(sprintf("Total transcripts: %d\n", nrow(utr_output)))
cat(sprintf("  5' UTR is contained within annotated boundaries (subset): %d\n", sum(utr_output$utr_cat == "subset", na.rm = TRUE)))
cat(sprintf("  5' UTR extends beyond annotated boundaries (unique): %d\n", sum(utr_output$utr_cat == "unique", na.rm = TRUE)))
cat(sprintf("  No open reading frame (no_orf): %d\n", sum(utr_output$utr_cat == "no_orf", na.rm = TRUE)))
cat(sprintf("  No 5' UTR (no_utr): %d\n", sum(utr_output$utr_cat == "no_utr", na.rm = TRUE)))


# =============================================================================
# STEP 4: SQANTI Protein classification
# =============================================================================

cat("\n--- STEP 4: SQANTI Protein classification ---\n")

protein_classifications = utr_output %>%
  mutate(
    protein_classification = case_when(
      
      ## Compared to the specific best-matching reference isoform
      # pr_nterm_diff: difference in ORF start positions (ATG)
      # pr_cterm_diff: difference in ORF end positions (stop codon)
      
      ## Compared to the nearest position across ALL isoforms of that gene
      # pr_nterm_gene_diff: distance to nearest CDS start among all isoforms
      # pr_cterm_gene_diff: distance to nearest CDS end among all isoforms
      
      # -----------------------------------------------------------------------------------------------------------
      # SQANTI_Protein, Full Splice Match - CDS splice junctions match exactly (known_splice), but precise N and C term may differ 
      # -----------------------------------------------------------------------------------------------------------
      
      # Case 1: Both ATG and stop codon match reference isoform exactly
      # FPM = Full protein match
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff == 0 & pr_cterm_diff == 0 ~ 
        "FPM,known_nterm_known_splice_known_cterm",
      
      # Case 2: Novel combination - ATG and/or STOP differ from reference isoform, but both exist in other isoforms of the gene
      # NPC = Novel protein combination
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        (pr_nterm_diff != 0 | pr_cterm_diff != 0) & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff == 0 ~ 
        "NPC,combo_nterm_cterm",
      
      # Case 3: ATG matches the reference exactly, but STOP codon is novel
      # NPE = Novel protein element
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff == 0 & pr_cterm_diff != 0 ~ 
        "NPE,known_nterm_known_splice_novel_cterm",
      
      # Case 4 subtypes: ATG doesn't match any isoform (Implicit: pr_nterm_gene_diff != 0), known STOP
      
      # 4a: Protein extends upstream - novel N-terminal amino acids
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nhang > 0 ~ 
        "NPE,novel_nterm_known_splice_known_cterm",
      
      # 4b: ATG is downstream/within reference CDS, using internal start codon = N-terminal truncation
      # IPM = Incomplete Protein Match
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nhang <= 0 & utr_cat == "subset" ~ 
        "IPM,nterm_truncation",
      
      # 4c: Novel TSS context - ATG position is novel but protein doesn't extend
      # Transcript has extended 5' UTR (novel TSS) creating novel regulatory context
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nhang <= 0 & utr_cat == "unique" ~ 
        "NPE,novel_nterm_known_splice_known_cterm",
      
      # Case 5 subtypes: Both ATG and stop codon differ from reference isoform
      
      # 5a: ATG matches another isoform (known), stop is novel
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff == 0 ~ 
        "NPE,known_nterm_known_splice_novel_cterm",
      
      # 5b: ATG doesn't match any isoform + internal start = truncation
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & utr_cat == "subset" ~ 
        "IPM,nterm_truncation",
      
      # 5c: Novel ATG (protein or TSS context) but end matches known isoform
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff == 0 ~ 
        "NPE,novel_nterm_known_splice_known_cterm",
      
      # 5d: Novel ATG and novel STOP
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff != 0 ~ 
        "NPE,novel_nterm_known_splice_novel_cterm",
      
      # FSM orphan catch-all
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" ~ 
        "orphan_fsm",
      
      # ---------------------------------------------------------------------------------------
      # SQANTI_Protein, Incomplete Splice Match (ISM) - CDS junctions are a subset of reference
      # (e.g., 5' fragment, 3' fragment, intron retention in CDS)
      # ---------------------------------------------------------------------------------------
      
      # Case 1: ATG matches reference isoform, stop codon is novel
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff == 0 & pr_cterm_diff != 0 & pr_cterm_gene_diff != 0 ~ 
        "NPE,known_nterm_known_splice_novel_cterm",
      
      # Case 2: Novel combination - ATG matches reference isoform, stop codon matches another isoform
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff == 0 & pr_cterm_diff != 0 & pr_cterm_gene_diff == 0 ~ 
        "NPC,combo_nterm_cterm",
      
      # Case 3: Novel combination - Both ATG and stop codon match (incomplete CDS splicing only)
      # Likely intron retention within CDS causing unique splicing pattern
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff == 0 & pr_cterm_diff == 0 ~ 
        "NPC,known_nterm_combo_splice_known_cterm",
      
      # Case 4: Novel combination - ATG differs from reference but matches another isoform, stop codon matches
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nterm_gene_diff == 0 ~ 
        "NPC,combo_nterm_cterm",
      
      # Case 5 subtypes - ATG doesnt match any isoform (novel), stop codon matches
      
      # Case 5a: ATG is downstream (internal start codon) - N-terminal truncation
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nterm_gene_diff != 0 & 
        pr_nhang <= 0 & utr_cat == "subset" ~ 
        "IPM,nterm_truncation",
      
      # Case 5b: Novel TSS context - ATG position novel but protein doesn't extend
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nterm_gene_diff != 0 & 
        pr_nhang <= 0 & utr_cat == "unique" ~ 
        "NPE,novel_nterm_known_splice_known_cterm",
      
      # Case 5c: Protein extends upstream - novel N-terminal amino acids
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nterm_gene_diff != 0 & pr_nhang > 0 ~ 
        "NPE,novel_nterm_known_splice_known_cterm",
      
      # Case 6: Novel Combination - Both ATG and stop match other isoforms
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff == 0 ~ 
        "NPC,combo_nterm_cterm",
      
      # Case 7: ATG matches another isoform, stop codon is novel
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff != 0 ~ 
        "NPE,known_nterm_known_splice_novel_cterm",
      
      # Case 8 subtypes - ATG doesn't match any isoform, stop codon varies
      
      # Case 8a: Internal start codon - truncation
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & utr_cat == "subset" ~ 
        "IPM,nterm_truncation",
      
      # Case 8b: Novel ATG (protein or TSS context), stop matches another isoform
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff == 0 ~ 
        "NPE,novel_nterm_known_splice_known_cterm",
      
      # Case 8c: Both ATG and stop codon are novel (don't match any isoform)
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff != 0 ~ 
        "NPE,novel_nterm_known_splice_novel_cterm",
      
      # ISM orphan catch-all
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") ~ 
        "orphan_ism",
      
      # ------------------------------------------------------------------------------------------
      # SQANTI_Protein, Novel in Catalog (NIC) - combo_splice (no reference isoform to compare to)
      # ------------------------------------------------------------------------------------------
      
      # Case 1: Both ATG and stop codon match an annotated isoform
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff == 0 & pr_cterm_gene_diff == 0 ~ 
        "NPC,known_nterm_combo_splice_known_cterm",
      
      # Case 2: ATG matches an isoform, stop codon is novel
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff == 0 & pr_cterm_gene_diff != 0 ~ 
        "NPE,known_nterm_combo_splice_novel_cterm",
      
      # Case 3: ATG doesn't match - truncation
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff != 0 & utr_cat == "subset" ~ 
        "IPM,nterm_truncation",
      
      # Case 4: Novel ATG (protein or TSS context), stop codon matches
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff != 0 & utr_cat == "unique" & pr_cterm_gene_diff == 0 ~ 
        "NPE,novel_nterm_combo_splice_known_cterm",
      
      # Case 5: Both ATG and stop codon are novel
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff != 0 & utr_cat == "unique" & pr_cterm_gene_diff != 0 ~ 
        "NPE,novel_nterm_combo_splice_novel_cterm",
      
      # NIC orphan catch-all
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") ~ 
        "orphan_nic",
      
      # ----------------------------------------------------------------------------------------------
      # SQANTI_Protein, Novel Not in Catalog (NNC) - novel_splice (no reference isoform to compare to)
      # ----------------------------------------------------------------------------------------------
      
      # Case 1: Both ATG and stop codon match an annotated isoform
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff == 0 ~ 
        "NPE,known_nterm_novel_splice_known_cterm",
      
      # Case 2: ATG matches an isoform, stop codon is novel
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff != 0 ~ 
        "NPE,known_nterm_novel_splice_novel_cterm",
      
      # Case 3: ATG doesn't match - truncation
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff != 0 & utr_cat == "subset" ~ 
        "IPM,nterm_truncation",
      
      # Case 4: Novel ATG (protein or TSS context), stop codon matches
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff == 0 ~ 
        "NPE,novel_nterm_novel_splice_known_cterm",
      
      # Case 5: Both ATG and stop codon are novel
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff != 0 ~ 
        "NPE,novel_nterm_novel_splice_novel_cterm",
      
      # NNC orphan catch-all
      pr_splice_cat == "novel_not_in_catalog" ~ 
        "orphan_nnc",
      
      # ------------------------------
      # Mono-exon cases
      # ------------------------------
      
      # Perfect mono-exon match
      (pr_splice_subcat == "mono-exon" | pr_splice_subcat == "mono-exon_by_intron_retention") & 
        pr_splice_cat == "full-splice_match" & pr_nterm_diff == 0 & pr_cterm_diff == 0 ~ 
        "FPM,mono-exon",
      
      # Intergenic mono-exon
      (pr_splice_subcat == "mono-exon" | pr_splice_subcat == "mono-exon_by_intron_retention") & 
        pr_splice_cat == "intergenic" ~ 
        "intergenic,mono-exon",
      
      # Other mono-exon cases (orphan)
      (pr_splice_subcat == "mono-exon" | pr_splice_subcat == "mono-exon_by_intron_retention") ~ 
        "orphan_monoexon,mono-exon",
      
      # ------------------------------
      # Miscellaneous multi-exon
      # ------------------------------
      
      pr_splice_subcat == "multi-exon" & pr_splice_cat == "intergenic" ~ 
        "intergenic,multi-exon",
      
      pr_splice_subcat == "multi-exon" & pr_splice_cat == "genic" ~ 
        "genic,multi-exon",
      
      pr_splice_subcat == "multi-exon" & pr_splice_cat == "antisense" ~ 
        "antisense,multi-exon",
      
      pr_splice_subcat == "multi-exon" & pr_splice_cat == "fusion" ~ 
        "fusion,multi-exon",
      
      # ------------------------------
      # Unclassified
      # ------------------------------
      
      TRUE ~ ""
    )
  ) %>%
  # Split protein_classification into base and subset
  separate(protein_classification, 
           into = c("protein_classification_base", "protein_classification_subset"), 
           sep = ",", 
           fill = "right", 
           remove = FALSE)

cat("\n--- Protein classification complete! ---\n")
cat("\nBy classification base:\n")
cat(sprintf("  Full protein match (FPM): %d\n", sum(protein_classifications$protein_classification_base == "FPM", na.rm = TRUE)))
cat(sprintf("  Incomplete protein match (IPM): %d\n", sum(protein_classifications$protein_classification_base == "IPM", na.rm = TRUE)))
cat(sprintf("  Novel protein combination (NPC): %d\n", sum(protein_classifications$protein_classification_base == "NPC", na.rm = TRUE)))
cat(sprintf("  Novel protein element (NPE): %d\n", sum(protein_classifications$protein_classification_base == "NPE", na.rm = TRUE)))
cat(sprintf("  intergenic: %d\n", sum(protein_classifications$protein_classification_base == "intergenic", na.rm = TRUE)))
cat(sprintf("  genic: %d\n", sum(protein_classifications$protein_classification_base == "genic", na.rm = TRUE)))
cat(sprintf("  antisense: %d\n", sum(protein_classifications$protein_classification_base == "antisense", na.rm = TRUE)))
cat(sprintf("  fusion: %d\n", sum(protein_classifications$protein_classification_base == "fusion", na.rm = TRUE)))
cat(sprintf("  orphan: %d\n", sum(grepl("orphan", protein_classifications$protein_classification_base), na.rm = TRUE)))

# =============================================================================
# STEP 5: Reconcile gene ids and gene names
# =============================================================================

cat("\n--- STEP 5: Reconcile gene ids and names ---\n")

gene_mapping = gencode_gtf %>% distinct(gene_id, gene_name)
ensg_lookup  = setNames(gene_mapping$gene_name, gene_mapping$gene_id)

# Just do direct assignment - no functions at all
protein_classifications$tx_gene_name = ensg_lookup[protein_classifications$tx_gene]
protein_classifications$pr_gene_name = ensg_lookup[protein_classifications$pr_gene]

# For the columns with commas, manually fix them (this is much faster than ifelse)
has_comma_tx = which(grepl(",", protein_classifications$tx_gene))
has_comma_pr = which(grepl(",", protein_classifications$pr_gene))

cat("Transcripts mapping to multiple genes at tx level: ", length(has_comma_tx), "\n")
cat("Transcripts mapping to multiple genes at pr level: ", length(has_comma_pr), "\n")

for (i in has_comma_tx) {
  ensgs = strsplit(protein_classifications$tx_gene[i], ",")[[1]]
  protein_classifications$tx_gene_name[i] = paste(unique(ensg_lookup[ensgs]), collapse = ",")
}

for (i in has_comma_pr) {
  ensgs = strsplit(protein_classifications$pr_gene[i], ",")[[1]]
  protein_classifications$pr_gene_name[i] = paste(unique(ensg_lookup[ensgs]), collapse = ",")
}

# Reconciliation of gene ids
matching = protein_classifications$tx_gene == protein_classifications$pr_gene
protein_classifications$reconciled_gene_id = protein_classifications$pr_gene  
protein_classifications$reconciled_gene_name = protein_classifications$pr_gene_name

# Only need to reconcile non-matching rows
non_matching = which(!matching)

cat(sprintf("\nSQANTI Transcript-level and protein-level gene assignments differ: %d of %d rows\n", 
            length(non_matching), nrow(protein_classifications)))

for (i in non_matching) {
  tx_ids  = strsplit(protein_classifications$tx_gene[i], ",")[[1]]
  pr_ids  = strsplit(protein_classifications$pr_gene[i], ",")[[1]]
  
  tx_name = strsplit(protein_classifications$tx_gene_name[i], ",")[[1]]
  pr_name = strsplit(protein_classifications$pr_gene_name[i], ",")[[1]]
  
  match_id   = intersect(tx_ids, pr_ids)
  match_name = intersect(tx_name, pr_name)
  
  if (length(match_id) > 0) {
    protein_classifications$reconciled_gene_id[i] = match_id[1]
  } else {
    protein_classifications$reconciled_gene_id[i] = pr_ids[1]
  }
  
  if (length(match_name) > 0) {
    protein_classifications$reconciled_gene_name[i] = match_name[1]
  } else {
    protein_classifications$reconciled_gene_name[i] = pr_name[1]
  }
}

# =============================================================================
# STEP 6: Filter transcripts post-SQANTI Protein
# =============================================================================

cat("\n--- STEP 6: Performing post-SQANTI Protein Filtering ---\n")

# Subclass lookup for short names, important for browser tracks
subclass_lookup = c(
  'known_nterm_novel_splice_known_cterm' = 'kn_ns_kc',
  'known_nterm_known_splice_known_cterm' = 'kn_ks_kc',
  'known_nterm_combo_splice_known_cterm' = 'kn_cs_kc',
  'known_nterm_novel_splice_novel_cterm' = 'kn_ns_nc',
  'known_nterm_combo_splice_novel_cterm' = 'kn_cs_nc',
  'known_nterm_known_splice_novel_cterm' = 'kn_ks_nc',
  'novel_nterm_known_splice_known_cterm' = 'nn_ks_kc',
  'novel_nterm_known_splice_novel_cterm' = 'nn_ks_nc',
  'novel_nterm_novel_splice_known_cterm' = 'nn_ns_kc',
  'novel_nterm_combo_splice_novel_cterm' = 'nn_cs_nc',
  'novel_nterm_combo_splice_known_cterm' = 'nn_cs_kc',
  'novel_nterm_novel_splice_novel_cterm' = 'nn_ns_nc',
  'combo_nterm_cterm'   = 'combo',
  'nterm_truncation'    = 'ntrunc',
  'multi-exon'          = 'multi',
  'mono-exon'           = 'mono'
)

full_protein = protein_classifications %>%
  left_join(exon_junctions, by = "transcript_id") %>%
  mutate(
    filter_status = case_when(
      
      # Check 1: Filter if no ORF found
      utr_cat == "no_orf" ~ "no_ORF",
      
      # Check 2: NMD filter - but rescue if only 1 junction and it's within 25bp of stop
      num_junc_after_stop_codon > min_junctions_after_stop_codon &
        !(num_junc_after_stop_codon == 1 & dist_stop_to_first_junc <= nmd_rescue_dist) ~ "NMD",
      
      # Check 3: Filter if problematic patterns
      grepl('intergenic|antisense|fusion|orphan|genic', protein_classification) ~ "sqanti_atypical",
      
      # Check 4: Filter if NOT in keep base classification list
      !(protein_classification_base %in% pclass_base_to_keep) ~ "sqanti_classification",
      
      # Keep everything else
      TRUE ~ "high_confidence"
    ),
    
    # abbreviate
    psubclass_short = subclass_lookup[protein_classification_subset]
    
  )

# NMD rescue
nmd_rescued = sum(
  full_protein$num_junc_after_stop_codon > min_junctions_after_stop_codon &
    full_protein$num_junc_after_stop_codon == 1 &
    full_protein$dist_stop_to_first_junc <= nmd_rescue_dist,
  na.rm = TRUE
)
cat(sprintf("NMD transcripts rescued (1 junc, <=%dbp from stop): %d\n", nmd_rescue_dist, nmd_rescued))

# Output full table with filter status
full_protein %<>%
  left_join(distinct(sample_gtf, gene_id, transcript_id), by = "transcript_id") %>%
  dplyr::select(transcript_id, gene_id, reference_gene_id = reconciled_gene_id, reference_gene_name = reconciled_gene_name, tx_transcripts, pr_transcripts, 
                num_junc_after_stop_codon, num_nt_after_stop_codon, num_5utr_exons, tss_in_gc_exons, utr_cat, 
                tclass = tx_cat, pclass = protein_classification_base, psubclass = protein_classification_subset, psubclass_short, filter_status)

write_tsv(full_protein, file.path(output_dir, paste0(basename, ".predicted_proteome.best_ORF_summary.txt")))

# Generate a matching amino acid fasta for mass spec reference- this will include lower confidence sequences such as NMD
full = readDNAStringSet(sample_dna_fasta_path) # dna fasta
full_fasta = tibble(
  transcript_id = names(full),
  full_dna_sequence = as.character(full)
)

# count matrix - include both counts and cpm columns for downstream analysis
counts = read_tsv(cpm_file_path) %>%
  select(transcript_id = 1, ends_with("_counts"), ends_with("_cpm"))

# orf coords
orf_info = read_tsv(mapped_orfs) %>%
  filter(orf_quality == "Clear Best ORF") %>%
  select(transcript_id = isoform_id, ORF_start, ORF_end) %>%
  left_join(counts, by = c("transcript_id")) %>%
  left_join(select(full_protein, 
                   gene_id, transcript_id, pclass, reference_gene_id, reference_gene_name, filter_status), 
            by = "transcript_id")

# group by orf
orf_groups = group_by_protein_sequence(orf_info, full_fasta)

# write fasta for all protein amino acid sequences for mass spec
orf_groups %<>%
  mutate(header = paste0(
    transcript_id,
    "|", gene_id, 
    "|", reference_gene_name,
    "|", pclass,
    "|", filter_status))

protein_seqs        = AAStringSet(orf_groups$orf_aa_sequence)
names(protein_seqs) = orf_groups$header
writeXStringSet(protein_seqs, file.path(output_dir, paste0(basename, ".predicted_proteome.best_ORF.fa")))

# =======================================================================================
# STEP 7: Write high confidence, ORF centric GTF and count file
# =======================================================================================

# table of high confidence orfs with cpm, for multisample analysis
hc_orf_info = orf_info %>% filter(filter_status == "high_confidence") 
hc_orf_groups = group_by_protein_sequence(hc_orf_info, full_fasta)

name_map = hc_orf_groups %>%
  group_by(orf_base_id) %>%
  summarise(reference_gene_name = paste(unique(reference_gene_name), collapse = ","), .groups = "drop")

hc_collapsed = hc_orf_groups %>%
  select(-avg_cpm) %>%
  group_by(orf_all_isoform_id, orf_base_id, gene_id) %>%
  summarize(across(c(ends_with("_counts"), ends_with("_cpm")), \(x) sum(x, na.rm = TRUE))) %>%
  ungroup() %>%
  left_join(name_map, by = "orf_base_id") %>%
  select(orf_base_id, orf_all_isoform_id, gene_id, reference_gene_name, everything())

write_tsv(hc_collapsed, file.path(output_dir, paste0(basename, ".predicted_proteome.collapsed_high_confidence_ORF_hashids_with_cpm.txt")))

# write a corresponding ORF centric gtf
gtf = import(sample_cds_gtf_path) %>%
  as.data.frame()

new_orf_attributes = hc_collapsed %>% 
  mutate(avg_orf_cpm = rowMeans(across(contains("cpm")), na.rm = TRUE)) %>%
  group_by(gene_id) %>%
  mutate(
    gene_total_cpm = sum(avg_orf_cpm, na.rm = TRUE),
    avg_orf_ratio = round(avg_orf_cpm / gene_total_cpm, 3)
  ) %>%
  ungroup() %>%
  select(orf_base_id, orf_all_isoform_id, reference_gene_name, avg_orf_ratio)

hc_gtf = gtf %>% 
  filter(type %in% c("transcript", "CDS")) %>%
  filter(transcript_id %in% hc_collapsed$orf_base_id)

orf_boundaries = hc_gtf %>%
  filter(type == "CDS") %>%
  group_by(transcript_id) %>%
  summarize(
    orf_start = min(start),
    orf_end = max(end)
  )

hc_gtf %<>% 
  left_join(orf_boundaries, by = "transcript_id") %>%
  mutate(
    start = if_else(type == "transcript" & !is.na(orf_start), orf_start, start),
    end   = if_else(type == "transcript" & !is.na(orf_end), orf_end, end)
  ) %>%
  select(-orf_start, -orf_end) %>%
  rename(orf_base_id = transcript_id) %>%
  left_join(new_orf_attributes, by = c("orf_base_id")) %>%
  mutate(name = paste0(orf_base_id, "|",
                       reference_gene_name, "|", 
                       avg_orf_ratio)) %>%
  select(-orf_all_isoform_id)

# calculate ratio
gr_updated = makeGRangesFromDataFrame(hc_gtf, keep.extra.columns = TRUE)
export(gr_updated, file.path(output_dir, paste0(basename, ".predicted_proteome.collapsed_high_confidence_ORF.gtf")), format = "gtf")

# convert to bed12
gtf_to_bed12(gtf_path = file.path(output_dir, paste0(basename, ".predicted_proteome.collapsed_high_confidence_ORF.gtf")),
             output_bed = file.path(output_dir, paste0(basename, ".predicted_proteome.collapsed_high_confidence_ORF.bed")),
             color_by = "avg_orf_ratio",
             track_name = paste0(basename, "_predicted_proteome"))
