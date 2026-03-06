#!/bin/bash

#SBATCH --job-name=protein_classification
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

BASENAME="all_AML"
DATA_DIR="/scratch/cwp5au/aml_test_data"
GENCODE_GTF="/project/sheynkman/external_data/GENCODE_v47/gencode.v47.annotation.gtf"

Rscript /scratch/cwp5au/LRP2_lite/bin/04_protein_classification.R \
  --basename ${BASENAME} \
  --gencode_gtf ${GENCODE_GTF} \
  --sample_cds_gtf ${DATA_DIR}/${BASENAME}_corrected_filtered_CDS.gtf \
  --sample_dna_fasta ${DATA_DIR}/${BASENAME}_corrected_filtered.fasta \
  --mapped_orfs ${DATA_DIR}/${BASENAME}_all_orfs_mapped.tsv \
  --protein_sqanti ${DATA_DIR}/${BASENAME}.sqanti_protein_classification.tsv \
  --cpm_file ${DATA_DIR}/${BASENAME}_hashids_with_cpm_filtered.txt \
  --output_dir ${DATA_DIR} \
  --min_junctions_after_stop 0 \
  --protein_class_keep "FPM,NPC,NPE" \
  --nmd_rescue_dist 25