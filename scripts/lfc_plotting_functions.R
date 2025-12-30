#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(magrittr)
library(patchwork)

options(scipen = 999)

# Updated by Megan Schertzer on December 18, 2025

# =============================================================================
# Function 1: Plot transcripts, allows multiple cluster ids
# =============================================================================

plot_transcripts <- function(exon_path, subisoforms_path, count_path, cluster_id = "1",
                             control_group = "WT", usage_threshold = 0.05) {
  
  subisoforms = fread(subisoforms_path)
  all_cluster_subisoforms = subisoforms %>% 
    filter(cluster_idx %in% cluster_id) %>%
    separate_rows(isoforms) 
  
  # Count total transcripts before filtering
  n_total_transcripts <- length(unique(all_cluster_subisoforms$isoforms))
  
  # Calculate transcript usage (proportion)
  transcript_usage = fread(count_path) %>%
    filter(isoform_id %in% all_cluster_subisoforms$isoforms) %>%
    select(gene_name, trans_id = isoform_id, ends_with("_cpm")) %>%
    group_by(gene_name) %>%
    mutate(across(ends_with("_cpm"), ~ . / sum(., na.rm = TRUE), .names = "{.col}_usage")) %>%
    ungroup() %>%
    rename_with(~ str_replace(., "_cpm_usage$", "_usage"), ends_with("_cpm_usage")) %>%
    mutate(
      max_usage = pmax(!!!select(., ends_with("_usage"))),
      avg_usage = rowMeans(select(., contains(control_group) & ends_with("_usage")), na.rm = TRUE)
    ) %>% 
    filter(max_usage >= usage_threshold) %>%
    arrange(avg_usage, trans_id) %>%
    select(gene_name, trans_id, avg_usage, max_usage)
  
  # Count transcripts after filtering
  n_removed <- n_total_transcripts - nrow(transcript_usage)
  
  # filter exons table
  exon_coords = fread(exon_path)
  gene_exon_coords = exon_coords %>% 
    filter(trans_id %in% transcript_usage$trans_id) %>%
    left_join(transcript_usage, by = "trans_id") %>%
    mutate(trans_id = factor(trans_id, levels = transcript_usage$trans_id))
  
  # for plotting
  for_plotting = gene_exon_coords %>% group_by(trans_id) %>% 
    summarise(start = min(start), end = max(end), .groups = "drop") 
  
  # shared axes info
  raw_limits   <- range(for_plotting$start, for_plotting$end)
  pad          <- diff(raw_limits)*0.05
  padded_limits<- raw_limits + c(-pad, +pad)
  
  shared_x <- scale_x_continuous(
    limits = padded_limits,
    expand = expansion(0, 0)
  )
  
  # plot transcripts
  p1 = ggplot() + 
    geom_segment(data = for_plotting, aes(x = start, xend = end, y = trans_id, yend = trans_id), size = 0.3) +
    geom_rect(data = gene_exon_coords, 
              aes(xmin = blockstarts, xmax = blockends, 
                  ymin = as.numeric(trans_id) -0.4, 
                  ymax = as.numeric(trans_id) + 0.4, 
                  fill = avg_usage)) +
    scale_fill_gradient(low = "lightsteelblue", high = "darkblue", name = "Control Usage") +
    geom_text(data = for_plotting, 
              aes(x = padded_limits[2], y = trans_id, label = trans_id), 
              hjust = 0, size = 3, check_overlap = TRUE) +
    shared_x + coord_cartesian(clip = "off") +
    theme_bw() +
    labs(title = sprintf("%s Transcript Summary (%s strand)", 
                         unique(gene_exon_coords$gene_name), 
                         unique(gene_exon_coords$strand)), 
         subtitle = sprintf("Showing %d/%d transcripts (removed %d with max usage < %.2f)",
                            nrow(transcript_usage), n_total_transcripts, n_removed, usage_threshold),
         x = "Genomic Position (bp)", y = NULL, fill = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5), 
          plot.margin = margin(0.1, 2, 0.1, 0.1, "in"), 
          legend.position = "left")
  
  return(list(plot = p1, limits = padded_limits, shared_x = shared_x))
}

# =============================================================================
# Function 2: Plot single cluster subisoforms
# =============================================================================

plot_cluster_subisoforms <- function(subisoforms_path, count_path, cluster_id = "1",
                                     control_group = "WT", usage_threshold = 0.05,
                                     padded_limits = NULL, shared_x = NULL) {
  
  subisoforms = fread(subisoforms_path)
  cluster_subisoforms = subisoforms %>% 
    filter(cluster_idx == cluster_id) %>%
    separate_rows(isoforms)
  
  # Calculate transcript usage (proportion) - same as plot_transcripts
  transcript_usage = fread(count_path) %>%
    filter(isoform_id %in% cluster_subisoforms$isoforms) %>%
    select(gene_name, trans_id = isoform_id, ends_with("_cpm")) %>%
    group_by(gene_name) %>%
    mutate(across(ends_with("_cpm"), ~ . / sum(., na.rm = TRUE), .names = "{.col}_usage")) %>%
    ungroup() %>%
    rename_with(~ str_replace(., "_cpm_usage$", "_usage"), ends_with("_cpm_usage")) %>%
    mutate(
      max_usage = pmax(!!!select(., ends_with("_usage"))),
      avg_usage = rowMeans(select(., contains(control_group) & ends_with("_usage")), na.rm = TRUE)
    ) %>% 
    filter(max_usage >= usage_threshold) %>%
    arrange(avg_usage, trans_id) %>%
    select(trans_id, avg_usage, max_usage)
  
  # Calculate subisoform usage considering ALL transcripts to match leafcutter output
  subisoform_coords = cluster_subisoforms %>% 
    select(-isoforms) %>%
    distinct() %>%
    group_by(cluster_idx) %>%
    mutate(across(ends_with("_counts"), ~ . / sum(., na.rm = TRUE), .names = "{.col}_proportion")) %>%
    ungroup() %>%
    mutate(avg_control = round(rowMeans(select(., contains(control_group) & ends_with("_proportion")), na.rm = TRUE), 2),
           avg_test = round(rowMeans(select(., !contains(control_group) & ends_with("_proportion")), na.rm = TRUE), 2))
  
  # Filter cluster_subisoforms to only include filtered transcripts
  cluster_subisoforms %<>%
    filter(isoforms %in% transcript_usage$trans_id)
  
  # Now filter subisoform_coords to match
  subisoform_coords %<>%
    filter(subisoform_id %in% cluster_subisoforms$subisoform_id)
  
  # Map subisoforms to isoforms
  map_to_isoform = cluster_subisoforms %>% 
    distinct(subisoform_id, isoforms) %>% 
    left_join(transcript_usage, by = c("isoforms" = "trans_id")) %>%
    arrange(desc(avg_usage), desc(isoforms)) %>%
    group_by(subisoform_id) %>%
    summarize(isoforms = paste0(isoforms, collapse = ",")) %>%
    ungroup()
  
  # Wrangle into blocks
  sub_introns = subisoform_coords %>%
    separate_rows(intron_starts, intron_ends, sep = ",", convert = TRUE) %>%
    group_by(subisoform_id) %>%
    mutate(
      junction_num = row_number(),
      n_junctions = n(),
      prev_intron_end = lag(intron_ends)
    ) %>%
    ungroup()
  
  plot_subs = bind_rows(
    
    # Middle exons
    sub_introns %>%
      filter(junction_num > 1) %>%
      mutate(block_start = prev_intron_end, block_end = intron_starts),
    
    # First anchor block
    sub_introns %>%
      filter(junction_num == 1) %>%
      mutate(block_start = intron_starts - 50, block_end = intron_starts),
    
    # Last anchor block
    sub_introns %>%
      filter(junction_num == n_junctions) %>%
      mutate(block_start = intron_ends, block_end = intron_ends + 50)
  ) %>%
    left_join(map_to_isoform, by = "subisoform_id") %>%
    mutate(plot_id = paste0("Control: ", avg_control, ", Test: ", avg_test))
  
  # Order by avg_control
  plot_id_order = plot_subs %>%
    distinct(isoforms, avg_control) %>%
    arrange(avg_control) %>%  
    pull(isoforms)
  
  plot_subs %<>%
    mutate(isoforms = factor(isoforms, levels = plot_id_order))
  
  # For plotting
  for_plotting_subs = plot_subs %>% 
    group_by(isoforms, plot_id) %>% 
    summarise(start = min(block_start), end = max(block_end), .groups = "drop")
  
  # Calculate limits and shared_x if not provided
  if (is.null(padded_limits)) {
    raw_limits <- range(for_plotting_subs$start, for_plotting_subs$end)
    pad        <- diff(raw_limits) * 0.05
    padded_limits <- raw_limits + c(-pad, +pad)
  }
  
  if (is.null(shared_x)) {
    shared_x <- scale_x_continuous(
      limits = padded_limits,
      expand = expansion(0, 0)
    )
  }
  
  # Plot subisoforms
  p = ggplot() + 
    geom_segment(data = for_plotting_subs, 
                 aes(x = start, xend = end, y = as.numeric(isoforms), 
                     yend = as.numeric(isoforms)), 
                 size = 0.3) +
    geom_rect(data = plot_subs, 
              aes(xmin = block_start, 
                  xmax = block_end, 
                  ymin = as.numeric(isoforms) -0.4, 
                  ymax = as.numeric(isoforms) + 0.4)) +
    geom_text(data = for_plotting_subs, 
              aes(y = as.numeric(isoforms), label = plot_id), 
              x = padded_limits[1], 
              hjust = 1, 
              size = 3, check_overlap = TRUE) +
    geom_text(data = for_plotting_subs, 
              aes(y = as.numeric(isoforms), label = isoforms), 
              x = padded_limits[2], 
              hjust = 0, 
              size = 3, check_overlap = TRUE) +
    shared_x + coord_cartesian(clip = "off") +
    theme_bw() +
    labs(title = paste("Cluster ", cluster_id), y = NULL, x = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          legend.position = "none", plot.margin = margin(0,1,0,0.1, "in")) 
  
  return(p)
}

# =============================================================================
# Function #3: Plot transcript usage
# =============================================================================

plot_transcript_usage <- function(subisoforms_path, count_path, cluster_id = "1", control_group = "WT", usage_threshold = 0.05) {
  
  # Get transcripts involved in clusters
  subisoforms = fread(subisoforms_path)
  all_cluster_subisoforms = subisoforms %>% 
    filter(cluster_idx %in% cluster_id) %>%
    separate_rows(isoforms) %>%
    distinct(isoforms)
  
  # Count total transcripts before filtering
  n_total_transcripts <- nrow(all_cluster_subisoforms)
  
  # Calculate transcript usage (proportion) - same as plot_transcripts
  transcript_usage = fread(count_path) %>%
    filter(isoform_id %in% all_cluster_subisoforms$isoforms) %>%
    select(gene_name, trans_id = isoform_id, ends_with("_cpm")) %>%
    group_by(gene_name) %>%
    mutate(across(ends_with("_cpm"), ~ . / sum(., na.rm = TRUE), .names = "{.col}_usage")) %>%
    ungroup() %>%
    rename_with(~ str_replace(., "_cpm_usage$", "_usage"), ends_with("_cpm_usage"))
  
  # Calculate gene-level CPM per sample (sum of transcript CPMs)
  gene_cpm = fread(count_path) %>%
    filter(isoform_id %in% all_cluster_subisoforms$isoforms) %>%
    select(trans_id = isoform_id, ends_with("_cpm")) %>%
    pivot_longer(cols = ends_with("_cpm"), names_to = "sample", values_to = "cpm") %>%
    mutate(sample = str_remove(sample, "_cpm")) %>%
    group_by(sample) %>%
    summarise(gene_cpm = sum(cpm, na.rm = TRUE), .groups = "drop") %>%
    mutate(group = ifelse(str_detect(sample, control_group), "Control", "Test"))
  
  # Calculate max and avg usage, filter and order
  trans_stats = transcript_usage %>%
    mutate(
      max_usage = pmax(!!!select(., ends_with("_usage"))),
      avg_usage = rowMeans(select(., contains(control_group) & ends_with("_usage")), na.rm = TRUE)
    ) %>%
    filter(max_usage >= usage_threshold) %>%
    arrange(desc(avg_usage), desc(trans_id)) %>%
    select(trans_id, max_usage, avg_usage)
  
  trans_order = trans_stats$trans_id
  
  # Count transcripts after filtering
  n_removed <- n_total_transcripts - length(trans_order)
  
  # Pivot to long format and filter
  transcript_usage_long = transcript_usage %>%
    filter(trans_id %in% trans_order) %>%
    select(gene_name, trans_id, ends_with("_usage")) %>%
    pivot_longer(cols = ends_with("_usage"), names_to = "sample", values_to = "usage") %>%
    mutate(sample = str_remove(sample, "_usage"),
           group = ifelse(str_detect(sample, control_group), "Control", "Test"),
           trans_id = factor(trans_id, levels = trans_order))
  
  # colors
  #n_trans <- length(unique(transcript_usage_long$trans_id))
  
  p = ggplot(transcript_usage_long, aes(x = sample, y = usage, fill = trans_id)) +
    geom_col(position = "stack") +
    geom_text(data = gene_cpm, 
              aes(x = sample, y = 1.05, label = paste0("CPM = ", round(gene_cpm, 0))), 
              inherit.aes = FALSE, size = 3, vjust = 0) +
    facet_wrap(~ group, scales = "free_x") +
    #scale_fill_gradientn(colors = colorRampPalette(c("lightsteelblue", "darkblue"))(n_trans), name = "") +
    scale_fill_viridis_d(name = "", option = "viridis") +
    #guides(fill = guide_legend(reverse = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw() +
    labs(title = paste(unique(transcript_usage$gene_name), "Transcript Usage"),
         subtitle = sprintf("Showing %d/%d transcripts (removed %d with max usage < %.2f)",
                 length(trans_order), n_total_transcripts, n_removed, usage_threshold),
         x = NULL, y = "Usage Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "top")
  
  return(p)
}

# =============================================================================
# Function #4: Plot subisoform usage
# =============================================================================

plot_subisoform_usage <- function(subisoforms_path, count_path, cluster_id = "1", control_group = "WT", usage_threshold = 0.05) {
  
  subisoforms = fread(subisoforms_path)
  cluster_subisoforms = subisoforms %>% 
    filter(cluster_idx == cluster_id) %>%
    separate_rows(isoforms)
  
  # Get list of filtered transcripts (same as plot_cluster_subisoforms)
  transcript_usage = fread(count_path) %>%
    filter(isoform_id %in% cluster_subisoforms$isoforms) %>%
    select(gene_name, trans_id = isoform_id, ends_with("_cpm")) %>%
    group_by(gene_name) %>%
    mutate(across(ends_with("_cpm"), ~ . / sum(., na.rm = TRUE), .names = "{.col}_usage")) %>%
    ungroup() %>%
    rename_with(~ str_replace(., "_cpm_usage$", "_usage"), ends_with("_cpm_usage")) %>%
    mutate(
      max_usage = pmax(!!!select(., ends_with("_usage"))),
      avg_usage = rowMeans(select(., contains(control_group) & ends_with("_usage")), na.rm = TRUE)
    ) %>% 
    filter(max_usage >= usage_threshold) %>%
    pull(trans_id)
  
  # Get subisoform usage across samples
  subisoform_usage_all = cluster_subisoforms %>%
    select(subisoform_id, ends_with("_counts")) %>%
    distinct() %>%
    pivot_longer(cols = ends_with("_counts"), names_to = "sample", values_to = "counts") %>%
    mutate(sample = str_remove(sample, "_counts")) %>%
    group_by(sample) %>%
    mutate(usage = counts / sum(counts, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(group = ifelse(str_detect(sample, control_group), "Control", "Test"))
  
  # Filter subisoforms based on transcript usage
  kept_subisoforms <- cluster_subisoforms %>%
    filter(isoforms %in% transcript_usage) %>%
    pull(subisoform_id) %>%
    unique()
  
  subisoform_usage <- subisoform_usage_all %>%
    filter(subisoform_id %in% kept_subisoforms)
  
  # Order subisoforms by average usage in control
  sub_order = subisoform_usage %>%
    filter(group == "Control") %>%
    group_by(subisoform_id) %>%
    summarise(avg_usage = mean(usage, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(avg_usage)) %>%
    pull(subisoform_id)
  
  subisoform_usage <- subisoform_usage %>%
    mutate(subisoform_id = factor(subisoform_id, levels = sub_order))
  
  p = ggplot(subisoform_usage, aes(x = sample, y = usage, fill = subisoform_id)) +
    geom_col(position = "stack") +
    facet_wrap(~ group, scales = "free_x") +
    scale_fill_viridis_d(name = "") +
    guides(fill = guide_legend(ncol = 1)) + 
    theme_bw() +
    labs(title = paste("Cluster", cluster_id, "Subisoform Usage"),
         x = NULL, y = "Usage Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "top")
  
  return(p)
}

# =============================================================================
# Main function to plot everything together
# =============================================================================

plot_transcript_subisoforms <- function(exon_path, subisoforms_path, count_path, cluster_id = "1", control_group = "WT", usage_threshold = 0.05, include_barplots = TRUE) {
  
  plots_list <- list()
  heights <- c()
  
  # 1. Transcript structure plot
  trans_result <- plot_transcripts(exon_path, subisoforms_path, count_path, cluster_id, control_group, usage_threshold)
  plots_list[[1]] <- trans_result$plot
  heights <- c(4)
  
  # 2. All cluster subisoform structure plots
  for (cid in cluster_id) {
    cluster_plot <- plot_cluster_subisoforms(subisoforms_path, count_path, cid, control_group, usage_threshold, 
                                             trans_result$limits, trans_result$shared_x)
    plots_list[[length(plots_list) + 1]] <- cluster_plot
    heights <- c(heights, 1)
  }
  
  # 3. Barplots layout depends on number of clusters
  if (include_barplots) {
    # Transcript usage barplot
    transcript_barplot <- plot_transcript_usage(subisoforms_path, count_path, cluster_id, control_group, usage_threshold)
    
    # Subisoform usage barplots
    subisoform_barplots <- list()
    for (cid in cluster_id) {
      usage_plot <- plot_subisoform_usage(subisoforms_path, count_path, cid, control_group, usage_threshold)
      subisoform_barplots[[length(subisoform_barplots) + 1]] <- usage_plot
    }
    
    n_clusters <- length(cluster_id)
    
    if (n_clusters <= 2) {
      # Side-by-side layout for 1-2 clusters
      subisoform_grid <- wrap_plots(subisoform_barplots, ncol = 2)
      barplots_row <- transcript_barplot | subisoform_grid
      plots_list[[length(plots_list) + 1]] <- barplots_row
      heights <- c(heights, 2)
    } else {
      # Stacked layout for 3+ clusters
      # Add transcript barplot
      plots_list[[length(plots_list) + 1]] <- transcript_barplot
      heights <- c(heights, 2)
      
      # Add subisoform grid below
      subisoform_grid <- wrap_plots(subisoform_barplots, ncol = 3)  # or ncol = 2
      plots_list[[length(plots_list) + 1]] <- subisoform_grid
      heights <- c(heights, 2)
    }
  }
  
  # Combine all plots vertically
  wrap_plots(plots_list, ncol = 1, heights = heights)
}

