#!/usr/bin/env Rscript

#' Protein Classification Script
#' 
#' Part 1: Classifies Long read derived 5' UTRs as either:
#'   - "subset" (N-terminal truncation, using internal start codon)
#'   - "unique" (novel start site, extending into new sequence)
#'   
#' Part 2: Protein Classification
#' Integrates SQANTI-protein output to classify CDS structures:
#'   - pFSM (Full Protein Match): CDS junctions exactly match reference
#'   - pISM (Incomplete Protein Match): CDS junctions are subset of reference
#'   - pNIC (Novel In Catalog): Novel combination of known N/splice/C
#'   - pNNC (Novel Not in Catalog): At least one novel N/splice/C

# =============================================================================
# Load libraries
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(txdbmaker)
  library(tidyverse)
  library(magrittr)
  library(rtracklayer)
})

# =============================================================================
# Get environment variables and check required files
# =============================================================================

basename           <- Sys.getenv("OUTPUT_BASE_NAME")
cpat_dir           <- file.path(Sys.getenv("OUTPUT_DIR"), "orf_calling") 
protein_sqanti_dir <- file.path(Sys.getenv("OUTPUT_DIR"), "protein_sqanti") 
gencode_gtf_path   <- Sys.getenv("GENCODE_GTF_FILE")

sample_cds_gtf_path   <- file.path(cpat_dir, paste0(basename, "_corrected_filtered_CDS.gtf"))
protein_sqanti_path   <- file.path(protein_sqanti_dir, paste0(basename, ".sqanti_protein_classification.tsv"))

stopifnot("GENCODE gtf file not found" = file.exists(gencode_gtf_path))
stopifnot("Sample gtf file not found" = file.exists(sample_cds_gtf_path))
stopifnot("Protein classification file not found" = file.exists(protein_sqanti_path))

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
  if (is.na(tss) || junc_chain == "" || length(gc_chains) == 0) return("novel")
  
  # Split the concatenated chains
  gc_chains = str_split(gc_chains, "\\|\\|\\|")[[1]]
  
  status = "novel"
  
  for (gc_chain in gc_chains) {
    if (!str_detect(gc_chain, fixed(junc_chain))) next
    
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

# =============================================================================
# STEP 1: Process GENCODE reference
# =============================================================================

message("\n--- STEP 1: Processing GENCODE reference ---")

# Load GENCODE GTF
gencode_gtf = import(gencode_gtf_path) %>% as.data.frame()

# Identify protein-coding transcripts (have CDS)
pc_transcripts = gencode_gtf %>%
  filter(type == "CDS") %>%
  pull(transcript_name) %>%
  unique()

# Extract exons and create GRanges for overlap checking
gencode_exons = gencode_gtf %>%
  filter(type == "exon", transcript_name %in% pc_transcripts) %>%
  select(seqnames, start, end, gene_id, transcript_name, strand)

gencode_gr = gencode_exons %>%
  select(seqnames, start, end, strand) %>%
  makeGRangesFromDataFrame()

gencode_gr = GenomicRanges::reduce(gencode_gr) # Reduces overlapping exons into merged genomic intervals

# Create exon chain strings per gene for junction matching (could potentially be tweaked to mimic collapse logic)
gencode_chains = gencode_exons %>%
  arrange(gene_id, transcript_name, start) %>%
  group_by(gene_id, transcript_name) %>%
  summarise(
    exon_chain = paste(paste0(start, "-", end), collapse = "_"),
    .groups = "drop_last"  # Keep gene grouping!
  ) %>%
  summarise(chains = paste(exon_chain, collapse = "|||"), .groups = "drop") 

# =============================================================================
# STEP 2: Process long read transcripts
# =============================================================================

message("\n--- STEP 2: Processing Sample transcripts ---")

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
    .groups = "drop"
  )

sample_exons = sample_gtf %>%
  filter(type == "exon") %>%
  group_by(transcript_id) %>%
  summarise(
    tss = ifelse(first(strand) == "+", min(start), max(end)),
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
# STEP 3: Classify 5' UTRs
# =============================================================================

message("\n--- STEP 3: Classifying 5' UTRs ---")

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
  left_join(gencode_chains, by = "gene_id")

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
  select(isoform_id = transcript_id, num_5utr_exons, utr_exon_status, tss_in_gc_exons, junc_cat, utr_cat)

# Merge with SQANTI protein classification
sqanti_class = read_tsv(protein_sqanti_path, show_col_types = FALSE)

# Adjust C term differences since CPAT includes stop codon in coords and gencode does not
sqanti_class %<>% mutate(
  pr_cterm_diff = ifelse(!is.na(pr_cterm_diff), pr_cterm_diff + 3, pr_cterm_diff),
  pr_cterm_gene_diff = ifelse(!is.na(pr_cterm_gene_diff), pr_cterm_gene_diff - 3, pr_cterm_gene_diff),
  pr_chang = ifelse(!is.na(pr_chang), pr_chang - 3, pr_chang))

utr_output = sqanti_class %>%
  full_join(utr_results, by = "isoform_id")

# write out intermediate file for troubleshooting
#write_tsv(utr_output, paste0("../results/protein_sqanti/", basename, "_5utr_classification.tsv"))


# 5' UTR summary statistics
message("\n--- 5'UTR Classification complete! ---")
message(sprintf("Total transcripts: %d", nrow(utr_output)))
message(sprintf("  subset: %d", sum(utr_output$utr_cat == "subset", na.rm = TRUE)))
message(sprintf("  unique: %d", sum(utr_output$utr_cat == "unique", na.rm = TRUE)))
message(sprintf("  no_orf: %d", sum(utr_output$utr_cat == "no_orf", na.rm = TRUE)))
message(sprintf("  no_utr: %d", sum(utr_output$utr_cat == "no_utr", na.rm = TRUE)))


# =============================================================================
# STEP 4: Protein classification
# =============================================================================

message("\n--- STEP 4: Protein classification ---")

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
        "pFPM,known_nterm_known_splice_known_cterm",
      
      # Case 2: Novel combination - ATG and/or STOP differ from reference isoform, but both exist in other isoforms of the gene
      # NPC = Novel protein combination
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        (pr_nterm_diff != 0 | pr_cterm_diff != 0) & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff == 0 ~ 
        "pNIC,combo_nterm_cterm",
      
      # Case 3: ATG matches the reference exactly, but STOP codon is novel
      # NPE = Novel protein element
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff == 0 & pr_cterm_diff != 0 ~ 
        "pNNC,known_nterm_known_splice_novel_cterm",
      
      # Case 4 subtypes: ATG doesn't match any isoform (Implicit: pr_nterm_gene_diff != 0), known STOP
      
      # 4a: Protein extends upstream - novel N-terminal amino acids
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nhang > 0 ~ 
        "pNNC,novel_nterm_known_splice_known_cterm",
      
      # 4b: ATG is downstream/within reference CDS, using internal start codon = N-terminal truncation
      # IPM = Incomplete Protein Match
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nhang <= 0 & utr_cat == "subset" ~ 
        "pIPM,nterm_truncation",
      
      # 4c: Novel TSS context - ATG position is novel but protein doesn't extend
      # Transcript has extended 5' UTR (novel TSS) creating novel regulatory context
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nhang <= 0 & utr_cat == "unique" ~ 
        "pNNC,novel_nterm_known_splice_known_cterm",
      
      # Case 5 subtypes: Both ATG and stop codon differ from reference isoform
      
      # 5a: ATG matches another isoform (known), stop is novel
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff == 0 ~ 
        "pNNC,known_nterm_known_splice_novel_cterm",
      
      # 5b: ATG doesn't match any isoform + internal start = truncation
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & utr_cat == "subset" ~ 
        "pIPM,nterm_truncation",
      
      # 5c: Novel ATG (protein or TSS context) but end matches known isoform
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff == 0 ~ 
        "pNNC,novel_nterm_known_splice_known_cterm",
      
      # 5d: Novel ATG and novel STOP
      pr_splice_cat == "full-splice_match" & pr_splice_subcat == "multi-exon" & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff != 0 ~ 
        "pNNC,novel_nterm_known_splice_novel_cterm",
      
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
        "pNNC,known_nterm_known_splice_novel_cterm",
      
      # Case 2: Novel combination - ATG matches reference isoform, stop codon matches another isoform
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff == 0 & pr_cterm_diff != 0 & pr_cterm_gene_diff == 0 ~ 
        "pNIC,combo_nterm_cterm",
      
      # Case 3: Novel combination - Both ATG and stop codon match (incomplete CDS splicing only)
      # Likely intron retention within CDS causing unique splicing pattern
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff == 0 & pr_cterm_diff == 0 ~ 
        "pNIC,known_nterm_combo_splice_known_cterm",
      
      # Case 4: Novel combination - ATG differs from reference but matches another isoform, stop codon matches
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nterm_gene_diff == 0 ~ 
        "pNIC,combo_nterm_cterm",
      
      # Case 5 subtypes - ATG doesnt match any isoform (novel), stop codon matches
      
      # Case 5a: ATG is downstream (internal start codon) - N-terminal truncation
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nterm_gene_diff != 0 & 
        pr_nhang <= 0 & utr_cat == "subset" ~ 
        "pIPM,nterm_truncation",
      
      # Case 5b: Novel TSS context - ATG position novel but protein doesn't extend
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nterm_gene_diff != 0 & 
        pr_nhang <= 0 & utr_cat == "unique" ~ 
        "pNNC,novel_nterm_known_splice_known_cterm",
      
      # Case 5c: Protein extends upstream - novel N-terminal amino acids
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff == 0 & pr_nterm_gene_diff != 0 & pr_nhang > 0 ~ 
        "pNNC,novel_nterm_known_splice_known_cterm",
      
      # Case 6: Novel Combination - Both ATG and stop match other isoforms
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff == 0 ~ 
        "pNIC,combo_nterm_cterm",
      
      # Case 7: ATG matches another isoform, stop codon is novel
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff != 0 ~ 
        "pNNC,known_nterm_known_splice_novel_cterm",
      
      # Case 8 subtypes - ATG doesn't match any isoform, stop codon varies
      
      # Case 8a: Internal start codon - truncation
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & utr_cat == "subset" ~ 
        "pIPM,nterm_truncation",
      
      # Case 8b: Novel ATG (protein or TSS context), stop matches another isoform
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff == 0 ~ 
        "pNNC,novel_nterm_known_splice_known_cterm",
      
      # Case 8c: Both ATG and stop codon are novel (don't match any isoform)
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_diff != 0 & pr_cterm_diff != 0 & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff != 0 ~ 
        "pNNC,novel_nterm_known_splice_novel_cterm",
      
      # ISM orphan catch-all
      pr_splice_cat == "incomplete-splice_match" & !str_detect(pr_splice_subcat, "mono-exon") ~ 
        "orphan_ism",
      
      # ------------------------------------------------------------------------------------------
      # SQANTI_Protein, Novel in Catalog (NIC) - combo_splice (no reference isoform to compare to)
      # ------------------------------------------------------------------------------------------
      
      # Case 1: Both ATG and stop codon match an annotated isoform
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff == 0 & pr_cterm_gene_diff == 0 ~ 
        "pNIC,known_nterm_combo_splice_known_cterm",
      
      # Case 2: ATG matches an isoform, stop codon is novel
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff == 0 & pr_cterm_gene_diff != 0 ~ 
        "pNNC,known_nterm_combo_splice_novel_cterm",
      
      # Case 3: ATG doesn't match - truncation
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff != 0 & utr_cat == "subset" ~ 
        "pIPM,nterm_truncation",
      
      # Case 4: Novel ATG (protein or TSS context), stop codon matches
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff != 0 & utr_cat == "unique" & pr_cterm_gene_diff == 0 ~ 
        "pNNC,novel_nterm_combo_splice_known_cterm",
      
      # Case 5: Both ATG and stop codon are novel
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") & 
        pr_nterm_gene_diff != 0 & utr_cat == "unique" & pr_cterm_gene_diff != 0 ~ 
        "pNNC,novel_nterm_combo_splice_novel_cterm",
      
      # NIC orphan catch-all
      pr_splice_cat == "novel_in_catalog" & !str_detect(pr_splice_subcat, "mono-exon") ~ 
        "orphan_nic",
      
      # ----------------------------------------------------------------------------------------------
      # SQANTI_Protein, Novel Not in Catalog (NNC) - novel_splice (no reference isoform to compare to)
      # ----------------------------------------------------------------------------------------------
      
      # Case 1: Both ATG and stop codon match an annotated isoform
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff == 0 ~ 
        "pNNC,known_nterm_novel_splice_known_cterm",
      
      # Case 2: ATG matches an isoform, stop codon is novel
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff == 0 & pr_cterm_gene_diff != 0 ~ 
        "pNNC,known_nterm_novel_splice_novel_cterm",
      
      # Case 3: ATG doesn't match - truncation
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff != 0 & utr_cat == "subset" ~ 
        "pIPM,nterm_truncation",
      
      # Case 4: Novel ATG (protein or TSS context), stop codon matches
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff == 0 ~ 
        "pNNC,novel_nterm_novel_splice_known_cterm",
      
      # Case 5: Both ATG and stop codon are novel
      pr_splice_cat == "novel_not_in_catalog" & pr_nterm_gene_diff != 0 & 
        utr_cat == "unique" & pr_cterm_gene_diff != 0 ~ 
        "pNNC,novel_nterm_novel_splice_novel_cterm",
      
      # NNC orphan catch-all
      pr_splice_cat == "novel_not_in_catalog" ~ 
        "orphan_nnc",
      
      # ------------------------------
      # Mono-exon cases
      # ------------------------------
      
      # Perfect mono-exon match
      (pr_splice_subcat == "mono-exon" | pr_splice_subcat == "mono-exon_by_intron_retention") & 
        pr_splice_cat == "full-splice_match" & pr_nterm_diff == 0 & pr_cterm_diff == 0 ~ 
        "pFPM,mono-exon",
      
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

message("\n--- Protein classification complete! ---")
message("\nBy classification base:")
message(sprintf("  pFPM: %d", sum(pb_classifications$protein_classification_base == "pFPM", na.rm = TRUE)))
message(sprintf("  pIPM: %d", sum(pb_classifications$protein_classification_base == "pIPM", na.rm = TRUE)))
message(sprintf("  pNIC: %d", sum(pb_classifications$protein_classification_base == "pNIC", na.rm = TRUE)))
message(sprintf("  pNNC: %d", sum(pb_classifications$protein_classification_base == "pNNC", na.rm = TRUE)))
message(sprintf("  intergenic: %d", sum(pb_classifications$protein_classification_base == "intergenic", na.rm = TRUE)))
message(sprintf("  genic: %d", sum(pb_classifications$protein_classification_base == "genic", na.rm = TRUE)))
message(sprintf("  antisense: %d", sum(pb_classifications$protein_classification_base == "antisense", na.rm = TRUE)))
message(sprintf("  fusion: %d", sum(pb_classifications$protein_classification_base == "fusion", na.rm = TRUE)))
message(sprintf("  orphan: %d", sum(grepl("orphan", pb_classifications$protein_classification_base), na.rm = TRUE)))

message("\nAnalysis complete!")