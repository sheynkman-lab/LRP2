#!/bin/bash
# config/config.sh - Pipeline Configuration File
# Edit this file with your specific paths and parameters

# =============================================================================
# General, relevant for entire pipeline
# =============================================================================

# Reference genome
export GENCODE_GENOME_FA="/home/cwp5au/datasets/jurkat_zenodo5234651/GRCh38.primary_assembly.genome.chr22.fa"
export GENCODE_GTF_FILE="/home/cwp5au/datasets/jurkat_zenodo5234651/gencode.v35.annotation.chr22.gtf"

# Output directory
export OUTPUT_DIR="/home/cwp5au/megan/github/LRP2_lite/jurkat_test_results"

# Base name for all output files in this project
export OUTPUT_BASE_NAME="jurkat_chr22"

# Species- currently supports human and mouse, but this is mainly a limitation of CPAT
export SPECIES="human"

# =============================================================================
# Environment
# =============================================================================

export LRP2_ENV_NAME="sqanti3"
export LRP2_ENV_FILE="environments/SQANTI3.conda_env.yml"

# =============================================================================
# Module: 01_Isoseq
# =============================================================================

# Input data paths
export SAMPLE_METADATA="config/jurkat_sample_metadata.txt"
export NUM_SAMPLES=1  # Number of *flnc.bam files in directory

# Iso-Seq collapse parameters- these are the default but the user can change
export MAX_FUZZY_JUNCTION=0
export MAX_5P_DIFF=100
export MAX_3P_DIFF=200

# =============================================================================
# Module: 02_sqanti
# =============================================================================

# SQANTI path- change this based on sqanti path
export SQANTI_PATH="/project/sheynkman/programs/SQANTI3-5.5"

# GTF source options: "lrp_isoseq" (this pipeline) or "custom" 
export GTF_SOURCE="lrp_isoseq"

# If using custom GTF, specify path here
export CUSTOM_GTF_PATH=""

# SQANTI filtering parameters
export PROTEIN_CODING_FILTER=TRUE
export INTERNAL_PRIMING_FILTER=TRUE
export TEMPLATE_SWITCHING_FILTER=TRUE
export PERCENT_POLYA_THRESHOLD=60
export TRANSCRIPT_CLASS_KEEP="FSM,NIC,NNC"

# =============================================================================
# Module: 03_orf_calling
# =============================================================================

# CPAT parameters, should we add --antisense?
export MIN_ORF=50 # cpat default is 75, same as NCBI ORFfinder
export TOP_ORF=50 # CPAT recommends 100

# =============================================================================
# Module: 04_sqanti_protein
# =============================================================================

export MIN_JUNCTIONS_AFTER_STOP_CODON=0 # parameter to identify NMD
export PROTEIN_CLASS_KEEP="FPM,NPC,NPE" # options are FPM,IPM,NPC,NPE, default only removes incomplete protein match (IPM)
<<<<<<< HEAD:conf/config_original/lrp2_jurkat_test_config.sh
=======





>>>>>>> main:config/lrp2_jurkat_test_config.sh
