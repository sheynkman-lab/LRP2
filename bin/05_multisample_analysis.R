#!/usr/bin/env Rscript

#' Multi-sample Analysis
#' 
#' Differential Gene, Transcript, and ORF expression with edgeR
#' Differential Usage at Transcript and ORF levels


# =============================================================================
# Load libraries
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(DRIMSeq)
  library(matrixStats)
  library(data.table)
  library(tidyverse)
})

options(scipen = 999)

# =============================================================================
# Get environment variables and check required files
# =============================================================================

basename                 <- Sys.getenv("OUTPUT_BASE_NAME")
sqanti_transcript_dir    <- file.path(Sys.getenv("OUTPUT_DIR"), "sqanti_transcript")
sqanti_protein_dir       <- file.path(Sys.getenv("OUTPUT_DIR"), "protein_sqanti")
multisample_analysis_dir <- file.path(Sys.getenv("OUTPUT_DIR"), "multisample_analysis")

count_file_path          <- file.path(sqanti_transcript_dir, paste0(basename, "_hashids_with_cpm_filtered.txt"))
orf_count_file_path      <- file.path(sqanti_protein_dir, paste0(basename, "_hashids_with_cpm_ORF.txt"))
metadata_file_path       <- Sys.getenv("SAMPLE_METADATA")

control_group      <- as.character(Sys.getenv("CONTROL_GROUP"))
experimental_group <- as.character(Sys.getenv("EXPERIMENTAL_GROUP"))

stopifnot("Transcript count file not found" = file.exists(count_file_path))
stopifnot("ORF count file not found" = file.exists(orf_count_file_path))
stopifnot("Sample sheet not found" = file.exists(metadata_file_path))

# =============================================================================
# Configure transcript-level data for edgeR
# =============================================================================

message("\n--- Reading transcript and ORF count files and sample sheet ---")

# Gene and transcript df
counts_raw         = as.data.frame(fread(count_file_path, header = TRUE))
count_cols         = grep("_counts$", colnames(counts_raw), value = TRUE)
cpm_cols           = grep("_cpm$", colnames(counts_raw), value = TRUE)
count_sample_names = sub("_counts$", "", count_cols)

# orf df- assumed that samples are in the same order as in the transcript level df above
orf_counts_raw         = as.data.frame(fread(orf_count_file_path, header = TRUE))
orf_count_cols         = grep("_counts$", colnames(orf_counts_raw), value = TRUE)
cpm_cols               = grep("_cpm$", colnames(orf_counts_raw), value = TRUE)
orf_count_sample_names = sub("_counts$", "", orf_count_cols)

# Check that transcript and ORF samples order matches
if (!identical(count_sample_names, orf_count_sample_names)) {
  stop("Sample order differs between transcript and ORF count files!\n",
       "Transcript samples: ", paste(count_sample_names, collapse = ", "), "\n",
       "ORF samples: ", paste(orf_count_sample_names, collapse = ", "))
}

# reorder sample sheet
sample_metadata = as.data.frame(fread(metadata_file_path, header = TRUE))

# Handle both 'name'/'group' and 'sample_name'/'condition' column naming conventions
if ("sample_name" %in% colnames(sample_metadata) && !("name" %in% colnames(sample_metadata))) {
  sample_metadata$name <- sample_metadata$sample_name
}
if ("condition" %in% colnames(sample_metadata) && !("group" %in% colnames(sample_metadata))) {
  sample_metadata$group <- sample_metadata$condition
}

sample_metadata = sample_metadata[match(count_sample_names, sample_metadata$name), ] # reorders

# Check for any mismatches
if (any(is.na(sample_metadata$name))) {
  stop("Some count columns don't have matching samples in metadata file!")
}

# When creating the factor, set control as reference (first level)
control   = control_group
experimental = experimental_group

group        = factor(sample_metadata$group, levels = c(control, experimental))
sample_names = sample_metadata$name

# Check for biological replicates?
#has_replicates = any(table(group) > 1)
#if (!has_replicates) {
#  cat("WARNING: No biological replicates detected. Using fixed BCV = 0.4 for differential analysis.\n")
#  cat("Results should be interpreted with caution as statistical power is very limited.\n\n")
#}

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

if (!identical(colnames(dge_raw$counts), dge_raw$samples$samples)) {
  stop("Count matrix columns and sample names don't match!")
}

# Filter, normalize, and estimate dispersion
dge = dge_raw[filterByExpr(dge_raw), , keep.lib.sizes = TRUE]
dge = calcNormFactors(dge, method = "TMM")
dge = estimateDisp(dge, robust = TRUE)

# Run differential expression test
result_dge  = exactTest(dge, pair = c(control, experimental)) # so that logFC>0 means higher in AML
dge_results = topTags(result_dge, n = Inf)$table

# Summary
cat(sprintf("Genes with FDR < 0.05: %d\n", sum(dge_results$FDR < 0.05)))
print(summary(decideTests(result_dge)))

# Save results- combine with gene ids and cpm results
gene_cpms = counts_raw %>%
  select(ensg_gene_id, gene_name, all_of(cpm_cols)) %>%
  group_by(ensg_gene_id, gene_name) %>%
  summarise(across(everything(), sum)) %>%
  ungroup() %>%
  distinct()

dge_results %>%
  rownames_to_column("ensg_gene_id") %>%
  left_join(gene_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DGE_edgeR_results.txt")))

cpm(dge_raw) %>%
  as.data.frame() %>%
  rownames_to_column("ensg_gene_id") %>%
  left_join(gene_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DGE_edgeR_raw_CPM_matrix.txt")))

cpm(dge) %>%
  as.data.frame() %>%
  rownames_to_column("ensg_gene_id") %>%
  left_join(gene_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DGE_edgeR_normalized_CPM_matrix.txt")))

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

if (!identical(colnames(dte_raw$counts), dte_raw$samples$samples)) {
  stop("Count matrix columns and sample names don't match!")
}

# Filter, normalize, and estimate dispersion
dte = dte_raw[filterByExpr(dte_raw), , keep.lib.sizes = TRUE]
dte = calcNormFactors(dte, method = "TMM")
dte = estimateDisp(dte, robust = TRUE)

# Run differential expression test
result_dte  = exactTest(dte, pair = c(control, experimental))
dte_results = topTags(result_dte, n = Inf)$table

# Summary
cat(sprintf("Transcripts with FDR < 0.05: %d\n", sum(dte_results$FDR < 0.05)))
print(summary(decideTests(result_dte)))

# Save results
transcript_cpms = counts_raw %>%
<<<<<<< HEAD:scripts/05_multisample_analysis.R
  select(isoform_id, hash_id, ensg_gene_id, gene_name, enst_transcript_id, transcript_name, all_of(cpm_cols)) %>%
=======
  select(isoform_id, hash_id, ensg_gene_id, gene_name, enst_transcript_id, any_of("transcript_name"), all_of(cpm_cols)) %>%
>>>>>>> dev:bin/05_multisample_analysis.R
  distinct()

dte_results %>%
  rownames_to_column("isoform_id") %>%
  left_join(transcript_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DTE_edgeR_results.txt")))

cpm(dte_raw) %>%
  as.data.frame() %>%
  rownames_to_column("isoform_id") %>%
  left_join(transcript_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DTE_edgeR_raw_CPM_matrix.txt")))

cpm(dte) %>%
  as.data.frame() %>%
  rownames_to_column("isoform_id") %>%
  left_join(transcript_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DTE_edgeR_normalized_CPM_matrix.txt")))

# Generate MD plot
pdf(file.path(multisample_analysis_dir, paste0(basename, "_DTE_MD_plot.pdf")), width = 8, height = 6)
plotMD(result_dte)
abline(h = c(-1, 1), col = "blue")
dev.off()

# =============================================================================
# Differential ORF Expression with edgeR
# =============================================================================
cat("Running ORF-level differential expression analysis...\n")

# Extract protein count matrix
counts_orf_matrix = orf_counts_raw %>%
  select(orf_isoform_id, all_of(orf_count_cols)) %>%
  column_to_rownames("orf_isoform_id") %>%
  as.data.frame()

# Create DGEList object
dpe_raw = DGEList(counts = counts_orf_matrix, group = group)

# Replace entire samples data frame from gene analysis, keep the original library sizes
dpe_raw$samples          = dge_raw$samples
colnames(dpe_raw$counts) = sample_names

if (!identical(colnames(dpe_raw$counts), dpe_raw$samples$samples)) {
  stop("Count matrix columns and sample names don't match!")
}

# Filter, normalize, and estimate dispersion
dpe = dpe_raw[filterByExpr(dpe_raw), , keep.lib.sizes = TRUE]
dpe = calcNormFactors(dpe, method = "TMM")
dpe = estimateDisp(dpe, robust = TRUE)

# Run differential expression test
result_dpe  = exactTest(dpe, pair = c(control, experimental))
dpe_results = topTags(result_dpe, n = Inf)$table

# Summary
cat(sprintf("ORFs with FDR < 0.05: %d\n", sum(dpe_results$FDR < 0.05)))
print(summary(decideTests(result_dpe)))

# Save results as TSV with IDs as columns
orf_cpms = orf_counts_raw %>%
  select(orf_isoform_id, hash_id, gene_id, gene_name, all_of(cpm_cols)) %>%
  distinct()

dpe_results %>%
  rownames_to_column("orf_isoform_id") %>%
  left_join(orf_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DE_ORF_edgeR_results.txt")))

cpm(dpe_raw) %>%
  as.data.frame() %>%
  rownames_to_column("orf_isoform_id") %>%
  left_join(orf_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DE_ORF_edgeR_raw_CPM_matrix.txt")))

cpm(dpe) %>%
  as.data.frame() %>%
  rownames_to_column("orf_isoform_id") %>%
  left_join(orf_cpms) %>%
  write_tsv(file.path(multisample_analysis_dir, paste0(basename, "_DE_ORF_edgeR_normalized_CPM_matrix.txt")))

# Generate MD plot
pdf(file.path(multisample_analysis_dir, paste0(basename, "_DE_ORF_MD_plot.pdf")), width = 8, height = 6)
plotMD(result_dpe)
abline(h = c(-1, 1), col = "blue")
dev.off()

# =============================================================================
# Differential Transcript Usage with DRIMSeq
# =============================================================================
cat("Running differential transcript usage analysis with DRIMSeq...\n")

counts_drimseq = counts_raw %>%
  select(gene_id = ensg_gene_id,      # DRIMSeq needs this named "gene_id"
         feature_id = isoform_id,     # DRIMSeq needs this named "feature_id"
         all_of(count_cols)) %>%
  rename_with(~ sub("_counts$", "", .), .cols = all_of(count_cols))

# Create sample table for DRIMSeq
sample_table = data.frame(
  sample_id = sample_names,
  group = group
)

# Create and filter DRIMSeq object
#counts_drimseq = head(x = as.data.frame(counts_drimseq), n = 200)
d = dmDSdata(counts = as.data.frame(counts_drimseq), samples = sample_table)

cat(sprintf("Before filtering: %d genes, %d transcripts\n", 
            length(unique(counts(d)$gene_id)), 
            nrow(counts(d))))

d = dmFilter(d, 
             min_samps_gene_expr    = nrow(sample_table)/2, # requires gene to be expressed in half of samples
             min_gene_expr          = 10, # based on gene level counts
             min_samps_feature_expr = 0, # isoform
             min_feature_expr       = 0) # isoform

cat(sprintf("After filtering: %d genes, %d transcripts\n", 
            length(unique(counts(d)$gene_id)), 
            nrow(counts(d))))

# Design, fit, and test
design_full = model.matrix(~ group, data = samples(d))
set.seed(123)
d = dmPrecision(d, design = design_full)
d = dmFit(d, design = design_full)
d = dmTest(d, coef = paste0("group", experimental))

# Extract results
dtu_gene_results = results(d) %>%
  dplyr::rename(lr_gene         = lr, 
                df_gene         = df,
                pvalue_gene     = pvalue, 
                adj_pvalue_gene = adj_pvalue)

dtu_tx_results = results(d, level = "feature") %>%
  dplyr::rename(isoform_id            = feature_id,
                lr_transcript         = lr,
                df_transcript        = df,
                pvalue_transcript     = pvalue,
                adj_pvalue_transcript = adj_pvalue)

# Get proportions from DRIMSeq and calculate means and deltas
transript_usage_table = transcript_cpms %>% 
  group_by(ensg_gene_id, gene_name) %>%
  mutate(across(ends_with("_cpm"), 
                ~ . / sum(., na.rm = TRUE),
                .names = "{.col}_prop")) %>%
  ungroup() %>%
  rename_with(~ str_replace(., "_cpm_prop$", "_prop"), ends_with("_cpm_prop"))

props = proportions(d) # these are calculated per group not per sample

props_long = props %>%
  pivot_longer(cols = -c(gene_id, feature_id), 
               names_to = "sample_id", 
               values_to = "proportion") %>%
  left_join(sample_table, by = "sample_id") %>% 
  distinct(gene_id, feature_id, proportion, group)

group_props = props_long %>%
  pivot_wider(names_from = group, 
              values_from = proportion,
              names_prefix = "group_prop_") %>%
  mutate(delta_proportion = .data[[paste0("group_prop_", levels(group)[2])]] - 
           .data[[paste0("group_prop_", levels(group)[1])]])

dtu_summary = dtu_tx_results %>%
  left_join(dtu_gene_results, by = "gene_id") %>%
  left_join(group_props, by = c("gene_id", "isoform_id" = "feature_id")) %>%
  left_join(transript_usage_table, by = c("gene_id" = "ensg_gene_id", "isoform_id")) %>%
<<<<<<< HEAD:scripts/05_multisample_analysis.R
  select(isoform_id, enst_transcript_id, transcript_name, hash_id, lr_transcript, pvalue_transcript, adj_pvalue_transcript, ends_with("_cpm"), ends_with("_prop"), starts_with("group_prop_"), delta_proportion,
=======
  select(isoform_id, enst_transcript_id, any_of("transcript_name"), hash_id, lr_transcript, pvalue_transcript, adj_pvalue_transcript, ends_with("_cpm"), ends_with("_prop"), starts_with("group_prop_"), delta_proportion,
>>>>>>> dev:bin/05_multisample_analysis.R
         gene_id, gene_name, lr_gene, df_gene, pvalue_gene, adj_pvalue_gene)
# Summary
cat(sprintf("Genes with DTU (adj_pvalue < 0.05): %d\n", 
            sum(dtu_gene_results$adj_pvalue_gene < 0.05, na.rm = TRUE)))
cat(sprintf("Transcripts with differential usage (adj_pvalue < 0.05): %d\n", 
            sum(dtu_tx_results$adj_pvalue_transcript < 0.05, na.rm = TRUE)))

# Save results
write_tsv(dtu_summary, file.path(multisample_analysis_dir, paste0(basename, "_DTU_transcript_DRIMSeq_summary.txt")))

# =============================================================================
# Differential ORF Usage with DRIMSeq
# =============================================================================
cat("Running differential ORF usage analysis with DRIMSeq...\n")

# Prepare DRIMSeq format from already-loaded ORF counts
counts_drimseq_orf = orf_counts_raw %>%
  select(gene_id,
         feature_id = orf_isoform_id,
         all_of(orf_count_cols)) %>%
  rename_with(~ sub("_counts$", "", .), .cols = all_of(orf_count_cols))

# Create sample table for DRIMSeq (reuse from transcript analysis)
sample_table_orf = data.frame(
  sample_id = sample_names,
  group = group
)

# Create and filter DRIMSeq object
d_orf = dmDSdata(counts = as.data.frame(counts_drimseq_orf), samples = sample_table_orf)

cat(sprintf("Before filtering (ORF DU): %d genes, %d ORFs\n", 
            length(unique(counts(d_orf)$gene_id)), 
            nrow(counts(d_orf))))

d_orf = dmFilter(d_orf, 
                 min_samps_gene_expr = nrow(sample_table_orf)/2,
                 min_gene_expr = 10,
                 min_samps_feature_expr = 0,
                 min_feature_expr = 0)

cat(sprintf("After filtering (ORF DU): %d genes, %d ORFs\n", 
            length(unique(counts(d_orf)$gene_id)), 
            nrow(counts(d_orf))))

# Design, fit, and test
design_full_orf = model.matrix(~ group, data = samples(d_orf))
set.seed(123)
d_orf = dmPrecision(d_orf, design = design_full_orf)
d_orf = dmFit(d_orf, design = design_full_orf)
d_orf = dmTest(d_orf, coef = paste0("group", experimental))

# Extract results
dpu_gene_results = results(d_orf) %>%
  dplyr::rename(lr_gene         = lr, 
                df_gene         = df,
                pvalue_gene     = pvalue, 
                adj_pvalue_gene = adj_pvalue)

dpu_orf_results = results(d_orf, level = "feature") %>%
  dplyr::rename(orf_isoform_id = feature_id,
                lr_orf         = lr,
                df_orf         = df,
                pvalue_orf     = pvalue,
                adj_pvalue_orf = adj_pvalue)

# Get proportions from DRIMSeq and calculate means and deltas
orf_usage_table = orf_cpms %>% 
  group_by(gene_id, gene_name) %>%
  mutate(across(ends_with("_cpm"), 
                ~ . / sum(., na.rm = TRUE),
                .names = "{.col}_prop")) %>%
  ungroup() %>%
  rename_with(~ str_replace(., "_cpm_prop$", "_prop"), ends_with("_cpm_prop"))

props_orf = proportions(d_orf)

props_orf_long = props_orf %>%
  pivot_longer(cols = -c(gene_id, feature_id), 
               names_to = "sample_id", 
               values_to = "proportion") %>%
  left_join(sample_table, by = "sample_id") %>% 
  distinct(gene_id, feature_id, proportion, group)

group_orf_props = props_orf_long %>%
  pivot_wider(names_from = group, 
              values_from = proportion,
              names_prefix = "group_prop_") %>%
  mutate(delta_proportion = .data[[paste0("group_prop_", levels(group)[2])]] - 
           .data[[paste0("group_prop_", levels(group)[1])]])

dpu_summary = dpu_orf_results %>%
  left_join(dpu_gene_results, by = "gene_id") %>%
  left_join(group_orf_props, by = c("gene_id", "orf_isoform_id" = "feature_id")) %>%
  left_join(orf_usage_table, by = c("gene_id", "orf_isoform_id")) %>%
  select(orf_isoform_id, hash_id, lr_orf, pvalue_orf, adj_pvalue_orf, ends_with("_cpm"), ends_with("_prop"), starts_with("group_prop_"), delta_proportion,
         gene_id, gene_name, lr_gene, df_gene, pvalue_gene, adj_pvalue_gene)

# Summary
cat(sprintf("Genes with DPU (adj_pvalue < 0.05): %d\n", 
            sum(dpu_gene_results$adj_pvalue_gene < 0.05, na.rm = TRUE)))
cat(sprintf("ORFs with differential usage (adj_pvalue < 0.05): %d\n", 
            sum(dpu_orf_results$adj_pvalue_orf < 0.05, na.rm = TRUE)))

# Save results
write_tsv(dpu_summary, file.path(multisample_analysis_dir, paste0(basename, "_DU_ORF_DRIMSeq_summary.txt")))