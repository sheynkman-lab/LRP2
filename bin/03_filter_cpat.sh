#!/bin/bash

#SBATCH --job-name=filter_CPAT
#SBATCH --time=02:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --mem=16G
#SBATCH --output=log_files/%x-%j.log
#SBATCH --error=log_files/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cwp5au@virginia.edu #your email address to receive notifications

module purge
module load gcc/11.4.0
module load openmpi/4.1.4
module load R/4.5.0

BASENAME=merged
DATA_DIR=/scratch/cwp5au/260306_lrp2_aml_test/
REFERENCE_GTF=/project/sheynkman/external_data/GENCODE_v47/gencode.v47.annotation.gtf
SQANTI_DIR=$DATA_DIR/results/S2_TRANSCRIPTOME/M3_FILTER_TRANSCRIPTOME
CPAT_DIR=$DATA_DIR/results/S3_PREDICTED_PROTEOME/M1_CPAT_ORF

Rscript /scratch/cwp5au/LRP2_lite/bin/03_filter_cpat.R \
  --basename $BASENAME \
  --cpat_fasta $CPAT_DIR/${BASENAME}_cpat.ORF_seqs.fa \
  --cpat_results $CPAT_DIR/${BASENAME}_cpat.ORF_prob.tsv \
  --sample_fasta $SQANTI_DIR/${BASENAME}.transcriptome.corrected_filtered.fasta \
  --sample_gtf $SQANTI_DIR/${BASENAME}.transcriptome.corrected_filtered.gtf \
  --reference_gtf $REFERENCE_GTF \
  --mapping_file $DATA_DIR/results/S2_TRANSCRIPTOME/M2_GENERATE_HASHIDS/${BASENAME}_hashids_mapping.txt \
  --output_dir $DATA_DIR/results/S3_PREDICTED_PROTEOME/M2_FILTER_CPAT