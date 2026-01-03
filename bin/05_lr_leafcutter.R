#!/usr/bin/env Rscript

#' Long read leafcutter
#' 
#' 

# =============================================================================
# Load libraries
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(magrittr)
  library(rtracklayer)
  library(Matrix) 
  library(igraph)
  library(optparse)
})

options(scipen = 999)

# =============================================================================
# Get environment variables and check required files
# =============================================================================

basename                 <- Sys.getenv("OUTPUT_BASE_NAME")
#sqanti_transcript_dir   <- file.path(Sys.getenv("OUTPUT_DIR"), "sqanti_transcript")
#sqanti_protein_dir      <- file.path(Sys.getenv("OUTPUT_DIR"), "protein_sqanti")
leafcutter_analysis_dir  <- file.path(Sys.getenv("OUTPUT_DIR"), "multisample_analysis", "lr_leafcutter")
metadata_file_path       <- Sys.getenv("SAMPLE_METADATA")

option_list <- list(
  make_option(c("--gtf"), type = "character", default = NULL, help = "GTF path"),
  make_option(c("--counts"), type = "character", default = NULL, help = "Count matrix path"),
  make_option(c("--mode"), type = "character", default = "exon", help = "CDS or exon")
)

opt <- parse_args(OptionParser(option_list = option_list))

transcript_gtf <- opt$gtf
count_file     <- opt$counts 
mode           <- opt$mode

stopifnot("Count file not found" = file.exists(count_file))
stopifnot("GTF file not found" = file.exists(transcript_gtf))
stopifnot("Sample sheet not found" = file.exists(metadata_file_path))


# =============================================================================
# Data wrangling helper scripts
# =============================================================================

# gtf to psl- originally adapted from FLAIR
gtf_to_psl = function(gtf_file_path, gtf_filter = "exon"){
  
  # Read the GTF file 
  message(paste0("Importing gtf and filtering by type = " , gtf_filter, "..."))
  gtf_data = import(gtf_file_path, format = "gtf") %>%
    as.data.table()
  
  exons_df = gtf_data %>% 
    filter(type == gtf_filter) %>%
    transmute(
      chrom = as.character(seqnames),
      start = as.numeric(start) - 1,
      end = as.numeric(end),
      strand = as.character(strand),
      transcript_id = as.character(transcript_id),
      gene_id = as.character(gene_id)
    )
  
  # Process by transcript
  message("Processing transcripts...")
  psl_data = exons_df %>%
    group_by(transcript_id, gene_id, chrom, strand) %>%
    summarise(
      starts = list(start),
      ends = list(end),
      .groups = "drop"
    ) %>%
    mutate(
      # Calculate block information
      blocksizes  = map2(starts, ends, ~ .y - .x),
      blockstarts = starts,
      
      # Reverse if needed (negative strand or descending coords)
      needs_reverse = map2_lgl(blockstarts, blocksizes, ~ length(.x) > 1 && .x[[1]] > .x[[2]]),
      blocksizes    = if_else(needs_reverse, map(blocksizes, rev), blocksizes),
      blockstarts   = if_else(needs_reverse, map(blockstarts, rev), blockstarts),
      
      # Calculate PSL fields
      blockcount = map_int(blockstarts, length),
      tstart     = map_dbl(blockstarts, ~ .x[[1]]),
      tend       = map2_dbl(blockstarts, blocksizes, ~ .x[[length(.x)]] + .y[[length(.y)]]),
      qsize      = map_dbl(blocksizes, ~ sum(unlist(.x))),
      qname      = paste0(transcript_id, "_", gene_id),
      
      # Calculate qstarts (cumulative positions)
      qstarts = map(blocksizes, ~ c(0, cumsum(head(unlist(.x), -1)))),
      
      # Convert lists to comma-separated strings
      blocksizes_str  = map_chr(blocksizes, ~ paste(c(unlist(.x), ""), collapse = ",")),
      qstarts_str     = map_chr(qstarts, ~ paste(c(.x, ""), collapse = ",")),
      blockstarts_str = map_chr(blockstarts, ~ paste(c(unlist(.x), ""), collapse = ","))
    )
  
  # Create final PSL format
  psl_output = psl_data %>%
    transmute(
      matches = 0, misMatches = 0, repMatches = 0, nCount = 0,
      qNumInsert = 0, qBaseInsert = 0, tNumInsert = 0, tBaseInsert = 0,
      strand = strand,
      qName  = qname,
      qSize  = qsize,
      qStart = 0,
      qEnd   = qsize,
      tName  = chrom,
      tSize  = 0,
      tStart = tstart,
      tEnd   = tend,
      blockCount = blockcount,
      blockSizes = blocksizes_str,
      qStarts = qstarts_str,
      tStarts = blockstarts_str
    )
  
  message("Done! Generated ", nrow(psl_output), " PSL entries")
  return(psl_output)

}

# psl to exon and intron centric tables- also from Isoviz
psl_to_coords <- function(psl_file_path){
  
  # input and format psl, fread is essential to read in the columns appropriately
  psl = fread(psl_file_path, header = FALSE, sep = "\t")
  psl %<>% distinct()
  
  df = psl %>% dplyr::select(chr = V14, start = V16, end = V17, id = V10, strand = V9, blocksizes = V19, blockstarts= V21) %>% 
    separate(col = id, into = c("trans_id", "gene_id"), sep = "_", extra = "merge", fill = "right") %>%
    mutate(blocksizes  = str_remove(blocksizes, ",$"), blockstarts = str_remove(blockstarts, ",$"))
  
  df_blocks = df %>% 
    separate_rows(blocksizes, blockstarts) # block coords are exons
  
  df_blocks$blockstarts = as.numeric(df_blocks$blockstarts)
  df_blocks$blocksizes  = as.numeric(df_blocks$blocksizes)
  df_blocks$blockends   = df_blocks$blockstarts + df_blocks$blocksizes
  
  # go from block coords to junction coords, split the data frame by trans_id
  list_of_data = split(df_blocks, df_blocks$trans_id)
  
  # Process each data frame in the list
  result = lapply(list_of_data, function(data){
    new_ends   = data$blockstarts[-1] + 1
    len        = nrow(data)
    new_starts = data$blockends[-len]
    n          = len - 1
    
    list(
      trans_id     = rep(data$trans_id[1], n),
      intron_start = new_starts,
      intron_end   = new_ends
    )
  })
  
  # Convert the list of lists back to three vectors
  trans_id     = unlist(lapply(result, `[[`, "trans_id"))
  intron_start = unlist(lapply(result, `[[`, "intron_start"))
  intron_end   = unlist(lapply(result, `[[`, "intron_end"))
  intron_data  = data.frame(trans_id, intron_start, intron_end)
  rownames(intron_data) = NULL
  
  trans_info = as.data.frame(df_blocks %>% dplyr::select(chr, trans_id, gene_id, strand) %>% 
                               distinct() %>% 
                               filter(chr != "chrY", chr != "chrM"))
  
  intron_data %<>% left_join(trans_info,  by = "trans_id") %>%
    dplyr::select(chr, intron_start, intron_end, gene_id, trans_id, strand)
  
  # create a coords id column
  intron_data %<>% unite(c("chr", "intron_start", "intron_end"), col = "junc_id", remove = FALSE, sep = "_")
  
  return(list(exon_coords = df_blocks, intron_coords = intron_data))
}

# get junction level counts from long read- need to be able to use this with multiple samples
isoform_to_junction <- function(intron_coords, counts){
  
  # reshape to long to handle multiple samples
  counts_long = counts %>%
    pivot_longer(
      cols = -trans_id,
      names_to = c("sample", ".value"),
      names_pattern = "(.+)_(counts|cpm)"
    ) %>%
    replace_na(list(counts = 0, cpm = 0))
  
  # join isoform->junction map, then aggregate to junction per sample
  jc = intron_coords %>% left_join(counts_long, by = "trans_id") %>% 
    group_by(junc_id, chr, intron_start, intron_end, gene_id, strand, sample) %>%
    dplyr::summarise(isoforms = paste0(unique(trans_id), collapse = ","), 
                     lr_junc_count = sum(counts, na.rm = TRUE),
                     lr_junc_cpm   = sum(cpm, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(chr, intron_start, junc_id, sample)
  
  return(jc)
  
}

# ==================================================================================
# Minicutter functions- adapted from Isoviz minicutter for clustering of long reads
# ==================================================================================

plot_cluster_sizes = function(juncs) {
  
  cluster_sizes = juncs %>% 
    group_by(cluster_idx) %>% 
    dplyr::summarize(n = n()) %>% 
    ungroup()
  
  ta = table(cluster_sizes$n)
  
  tibble(
    cluster_size = as.numeric(names(ta)),
    num_clusters = as.numeric(ta)) %>%
    ggplot(aes(cluster_size, num_clusters)) + geom_point() + xlab("cluster size") + ylab("number clusters") + scale_y_log10()
}

calculate_usage_ratios = function(juncs) {
  juncs %>% group_by(cluster_idx) %>%
    mutate(usage_ratio = readcount / sum(readcount)) %>%
    ungroup() # calculate usage ratios
}

leafcutter_one_step = function(juncs) {
  
  juncs = juncs %>% 
    dplyr::select("chrom", "strand", "start", "end", "name", "readcount")
  
  splice_sites = bind_rows(
    juncs %>% dplyr::select(chrom, strand, position = start),
    juncs %>% dplyr::select(chrom, strand, position = end)) %>%
    distinct() %>%
    arrange(chrom, strand, position) %>%
    dplyr::mutate(idx = 1:n())
  
  juncs = juncs %>%
    left_join(splice_sites, by = c(chrom = "chrom", strand = "strand", start = "position")) %>%
    left_join(splice_sites, by = c(chrom = "chrom", strand = "strand", end = "position"),
              suffix = c("_start","_end"))
  
  nss = nrow(splice_sites)
  
  intron_connectivity <- sparseMatrix(
    i         = juncs$idx_start,
    j         = juncs$idx_end,
    dims      = c(nss, nss),
    x         = 1,
    symmetric = TRUE
  )
  
  g = graph_from_adjacency_matrix(intron_connectivity, "undirected")
  juncs$cluster_idx = igraph::components(g)$membership[juncs$idx_start]
  return(juncs)
}

# main function- default parameters keep almost all junctions. Can increase min_usage_ratio to remove low count junctions.
minicutter = function(juncs, plot_summary = TRUE, min_usage_ratio = 0.01) {
  
  colnames(juncs) = c("chrom", "start", "end", "name", "readcount", "strand")
  juncs           = leafcutter_one_step(juncs)
  
  if(plot_summary){
    print("Printing summary of intron-cluster sizes = how many junctions across in clusters")
    p = plot_cluster_sizes(juncs)
    print(p)
  }
  
  juncs = juncs %>% 
    calculate_usage_ratios()
  
  # refine clusters 
  juncs_filtered  = as.data.table(juncs %>% filter(usage_ratio >= min_usage_ratio))
  juncs_recluster = leafcutter_one_step(juncs_filtered) 
  
  if(plot_summary){
    print("Printing summary of intron-cluster sizes = how many junctions across in clusters")
    p = plot_cluster_sizes(juncs_filtered) + ggtitle("Post cluster refinemenet and removing lowly used junctions")
    print(p)
  }
  
  return(juncs_recluster)
}

# ==================================================================================
# Main
# ==================================================================================

# Create and output psl
psl_base = tools::file_path_sans_ext(basename(transcript_gtf))
psl_path = file.path(leafcutter_analysis_dir, paste0(psl_base, ".psl"))
gtf_to_psl(transcript_gtf, gtf_filter = mode) %>%
  write_tsv(psl_path, col_names = FALSE)

coords        = psl_to_coords(psl_path)
intron_coords = coords$intron_coords
exon_coords   = coords$exon_coords

write_tsv(intron_coords, file.path(leafcutter_analysis_dir, paste0(psl_base, "_intron_coords.txt")))
write_tsv(exon_coords, file.path(leafcutter_analysis_dir, paste0(psl_base, "_exon_coords.txt")))

# Collapse long read counts to junctions
counts = fread(count_file) %>% 
  select(trans_id = 1, ends_with("_counts"), ends_with("_cpm"))

lr_junctions = isoform_to_junction(intron_coords, counts)

# get single counts per junction- sum across samples and get in format for minicutter
lr_junctions_collapsed = lr_junctions %>% 
  group_by(chr, intron_start, intron_end, junc_id, strand) %>%
  summarize(readcount = sum(lr_junc_count, na.rm = TRUE)) %>%
  ungroup() %>%
  relocate(strand, .after = last_col())

# Run minicutter and select non-constitutive junctions
clusters = minicutter(juncs = lr_junctions_collapsed, plot_summary = TRUE, min_usage_ratio = 0.01)
as_clusters = clusters %>% 
  select(chr = chrom, start, end, junc_id = name, cluster_idx, strand) %>%
  group_by(cluster_idx) %>%
  filter(n() > 1) %>%
  ungroup()

# Overlap junctions/clusters with isoforms
junc_to_iso = lr_junctions %>% distinct(junc_id, isoforms)
by_isoform = as_clusters %>%
  left_join(junc_to_iso, by = c("junc_id")) %>%
  separate_rows(isoforms, sep = ",") %>%
  group_by(chr, cluster_idx, isoforms, strand) %>%
  arrange(start, end) %>%
  summarise(intron_starts = paste0(unique(start), collapse = ","),
            intron_ends = paste0(unique(end), collapse = ",")) %>%
  ungroup() 

# Need to sum per sample per subisoform
by_isoform %<>% 
  left_join(select(counts, -ends_with("_cpm")), by = c("isoforms" = "trans_id"))

subisoform = by_isoform %>%
  group_by(chr, cluster_idx, strand, intron_starts, intron_ends) %>%
  summarise(isoforms = paste0(unique(isoforms), collapse = ","),
            across(ends_with("_counts"), sum),
            .groups = "drop") %>%
  mutate(subisoform_id = paste0(chr, ":", intron_starts, ":", intron_ends, ":clu_", cluster_idx)) %>%
  relocate((subisoform_id))

write_tsv(subisoform, file = file.path(leafcutter_analysis_dir, paste0(basename, "_lr_leafcutter_subisoform_clusters.txt")), col_names = TRUE)

# write out space delimited file with row names for diff splicing script
out = subisoform %>%
  select(subisoform_id, ends_with("_counts")) %>%
  rename_with(~gsub("_counts$", "", .), ends_with("_counts")) %>%
  column_to_rownames("subisoform_id")

write.table(out, 
            file = gzfile(file.path(leafcutter_analysis_dir, paste0(basename, "_lr_leafcutter_perind_numers.counts.gz"))),
            sep = " ",
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE)

# create the group text file required for ds (use sample sheet)
sample_metadata = as.data.frame(fread(metadata_file_path, header = TRUE))

# Handle both 'name'/'group' and 'sample_name'/'condition' column naming conventions
if ("sample_name" %in% colnames(sample_metadata) && !("name" %in% colnames(sample_metadata))) {
  sample_metadata$name <- sample_metadata$sample_name
}
if ("condition" %in% colnames(sample_metadata) && !("group" %in% colnames(sample_metadata))) {
  sample_metadata$group <- sample_metadata$condition
}

sample_metadata %>% select(name, group) %>%
  write_tsv(file.path(leafcutter_analysis_dir, paste0(basename, "_groups_file.txt")), col_names = FALSE)