#!/usr/bin/env Rscript

#' Multi-sample Analysis
#' 
#' Differential Gene expression with edgeR


# =============================================================================
# Load libraries
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
  library(tidyverse)
})

# =============================================================================
# Get environment variables and check required files
# =============================================================================

basename                        <- Sys.getenv("OUTPUT_BASE_NAME")
sqanti_transcript_dir           <- file.path(Sys.getenv("OUTPUT_DIR"), "sqanti_transcript")
multisample_analysis_dir        <- file.path(Sys.getenv("OUTPUT_DIR"), "multisample_analysis")

count_file_path         <- file.path(sqanti_transcript_dir, paste0(basename, "_hashids_with_cpm_filtered.txt"))
metadata_file_path      <- Sys.getenv("SAMPLE_METADATA")

stopifnot("Count file not found" = file.exists(count_file_path))
stopifnot("Sample sheet not found" = file.exists(metadata_file_path))

# =============================================================================
# Configure data for edgeR
# =============================================================================

message("\n--- Reading count file and sample sheet ---")

counts_raw      = as.data.frame(fread(count_file_path, header = TRUE))
sample_metadata = as.data.frame(fread(metadata_file_path, header = TRUE))

# Extract count columns (those ending with "_counts")
count_cols = grep("_counts$", colnames(counts_raw), value = TRUE)
count_sample_names = sub("_counts$", "", count_cols)
sample_metadata = sample_metadata[match(count_sample_names, sample_metadata$name), ] # reorders

# Check for any mismatches
if (any(is.na(sample_metadata$name))) {
  stop("Some count columns don't have matching samples in metadata file!")
}

group        = factor(sample_metadata$group)
sample_names = sample_metadata$name

# =============================================================================
# Differential Gene Expression with edgeR
# =============================================================================

cat("Running gene-level differential expression analysis...\n")

# Aggregate transcripts to genes using ensg_gene_id
counts_gene_matrix = counts_raw %>%
  select(ensg_gene_id, all_of(count_cols)) %>%
  group_by(ensg_gene_id) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("ensg_gene_id") %>%
  as.data.frame()

# Create DGEList object
dge_raw = DGEList(counts = counts_gene_matrix, group = group)
colnames(dge_raw$counts) = sample_names
dge_raw$samples$samples  = sample_names

# Filter, normalize, and estimate dispersion
dge = dge_raw[filterByExpr(dge_raw), , keep.lib.sizes = FALSE]
dge = calcNormFactors(dge, method = "TMM")
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)

# Run differential expression test
result_dge  = exactTest(dge)
dge_results = topTags(result_dge, n = Inf)$table

# Summary
cat(sprintf("Genes with FDR < 0.05: %d\n", sum(dge_results$FDR < 0.05)))
print(summary(decideTests(result_dge)))

# Save results- combine with gene ids and cpm results
cpm_cols = grep("_cpm$", colnames(counts_raw), value = TRUE)
gene_cpms = counts_raw %>%
  select(ensg_gene_id, gene_name, all_of(cpm_cols)) %>%
  group_by(ensg_gene_id, gene_name) %>%
  summarise(across(everything(), sum)) %>%
  ungroup() %>%
  distinct()

dge_results %>%
  rownames_to_column("ensg_gene_id") %>%
  left_join(gene_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DGE_results.txt")))

cpm(dge_raw) %>%
  as.data.frame() %>%
  rownames_to_column("ensg_gene_id") %>%
  left_join(gene_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DGE_raw_CPM_matrix.txt")))

cpm(dge) %>%
  as.data.frame() %>%
  rownames_to_column("ensg_gene_id") %>%
  left_join(gene_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DGE_normalized_CPM_matrix.txt")))

# Generate MD plot
pdf(file.path(multisample_analysis_dir, paste0(basename, "_DGE_MD_plot.pdf")), width = 8, height = 6)
plotMD(result_dge)
abline(h = c(-1, 1), col = "blue")
dev.off()

# =============================================================================
# Differential Transcript Expression with edgeR
# =============================================================================
cat("Running transcript-level differential expression analysis...\n")

# Extract transcript count matrix
counts_transcript_matrix = counts_raw %>%
  select(isoform_id, all_of(count_cols)) %>%
  column_to_rownames("isoform_id") %>%
  as.data.frame()

# Create DGEList object
dte_raw = DGEList(counts = counts_transcript_matrix, group = group)
colnames(dte_raw$counts) = sample_names
dte_raw$samples$samples  = sample_names

# Filter, normalize, and estimate dispersion
dte = dte_raw[filterByExpr(dte_raw), , keep.lib.sizes = FALSE]
dte = calcNormFactors(dte, method = "TMM")
dte = estimateCommonDisp(dte)
dte = estimateTagwiseDisp(dte)

# Run differential expression test
result_dte  = exactTest(dte)
dte_results = topTags(result_dte, n = Inf)$table

# Summary
cat(sprintf("Transcripts with FDR < 0.05: %d\n", sum(dte_results$FDR < 0.05)))
print(summary(decideTests(result_dte)))

# Save results
transcript_cpms = counts_raw %>%
  select(isoform_id, hash_id, ensg_gene_id, gene_name, enst_transcript_id, transcript_name, all_of(cpm_cols)) %>%
  distinct()

dte_results %>%
  rownames_to_column("isoform_id") %>%
  left_join(transcript_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DTE_results.txt")))

cpm(dte_raw) %>%
  as.data.frame() %>%
  rownames_to_column("isoform_id") %>%
  left_join(transcript_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DTE_raw_CPM_matrix.txt")))

cpm(dte) %>%
  as.data.frame() %>%
  rownames_to_column("isoform_id") %>%
  left_join(transcript_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DTE_normalized_CPM_matrix.txt")))

# Generate MD plot
pdf(file.path(multisample_analysis_dir, paste0(basename, "_DTE_MD_plot.pdf")), width = 8, height = 6)
plotMD(result_dte)
abline(h = c(-1, 1), col = "blue")
dev.off()
