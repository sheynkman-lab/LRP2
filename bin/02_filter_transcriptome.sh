#!/bin/bash

#SBATCH --job-name=filter_transcriptome
#SBATCH --time=02:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --mem=32G
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
SQANTI_DIR=$DATA_DIR/results/S2_TRANSCRIPTOME/M1_SQANTI_QC

Rscript /scratch/cwp5au/LRP2_lite/bin/02_filter_sqanti_transcripts.R \
  --basename $BASENAME \
  --classification $SQANTI_DIR/${BASENAME}_classification.txt \
  --sample_gtf $SQANTI_DIR/${BASENAME}_corrected.gtf \
  --sample_fasta $SQANTI_DIR/${BASENAME}_corrected.fasta \
  --mapping_file $DATA_DIR/results/S2_TRANSCRIPTOME/M2_GENERATE_HASHIDS/${BASENAME}_hashids_mapping.txt \
  --output_dir $DATA_DIR/results/S2_TRANSCRIPTOME/M3_FILTER_TRANSCRIPTOME
  