#!/bin/bash
# config/config.sh - Pipeline Configuration File
# Edit this file with your specific paths and parameters

# =============================================================================
# Cluster info
# =============================================================================

export ACCOUNT="your_account_name"
export EMAIL="your_email@university.edu"

# =============================================================================
# User-specific paths
# =============================================================================

# Input data paths
export SAMPLE_METADATA="config/sample_metadata.txt"
export NUM_SAMPLES=4  # Number of *flnc.bam files in directory

# Reference genome
export GENOME_FA="/project/sheynkman/external_data/GENCODE_v47/GRCh38.primary_assembly.genome.fa"
export GTF_FILE="/project/sheynkman/external_data/GENCODE_v47/gencode.v47.annotation.gtf"

# Output directory
export OUTPUT_DIR="/home/cwp5au/megan/github/LRP2_lite/results"

# =============================================================================
# Tool parameters
# =============================================================================

# Iso-Seq collapse parameters
export MAX_FUZZY_JUNCTION=0
export MAX_5P_DIFF=100
export MAX_3P_DIFF=200