#!/bin/bash

#SBATCH --job-name=generate_hashids
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
REFERENCE_GTF="/project/sheynkman/external_data/GENCODE_v47/gencode.v47.annotation.gtf"

Rscript /scratch/cwp5au/LRP2_lite/bin/00_generate_hashids.R \
  --basename $BASENAME \
  --sample_gtf $DATA_DIR/results/S2_TRANSCRIPTOME/M1_SQANTI_QC/${BASENAME}_corrected.gtf \
  --classification $DATA_DIR/results/S2_TRANSCRIPTOME/M1_SQANTI_QC/${BASENAME}_classification.txt \
  --reference_gtf $REFERENCE_GTF \
  --hashlib_script /scratch/cwp5au/LRP2_lite/bin/00_hashlib_id_generator.py \
  --output_dir $DATA_DIR/results/S2_TRANSCRIPTOME/M2_GENERATE_HASHIDS
  