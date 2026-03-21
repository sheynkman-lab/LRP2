#!/bin/bash

#SBATCH --job-name=multisample
#SBATCH --time=02:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=log_files/%x-%j.log
#SBATCH --error=log_files/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cwp5au@virginia.edu #your email address to receive notifications

module purge
module load gcc/11.4.0
module load openmpi/4.1.4
module load R/4.5.0

BASENAME="merged"
DATA_DIR="/scratch/cwp5au/260306_lrp2_aml_test/"

Rscript /scratch/cwp5au/LRP2_lite/bin/05_multisample_analysis.R \
  --basename $BASENAME \
  --count_file ${DATA_DIR}/results/S2_TRANSCRIPTOME/M3_FILTER_TRANSCRIPTOME/${BASENAME}_transcriptome_hashids_with_cpm_filtered.txt \
  --orf_count_file ${DATA_DIR}/results/S3_PREDICTED_PROTEOME/M4_PROTEIN_UTR_CLASSIFICATION/${BASENAME}.predicted.proteome.high_confidence_ORF_cpm.txt \
  --sample_metadata ${DATA_DIR}/samplesheet_AML.csv \
  --control_group "NBM" \
  --experimental_group "AML" \
  --output_dir ${DATA_DIR}/results/S4_MULTISAMPLE_ANALYSIS \
  --threads 8