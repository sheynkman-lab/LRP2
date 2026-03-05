#!/bin/bash

#SBATCH --job-name=build_reference
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

# Specify GENCODE version and output directory - gencode downloaded to output directory
# this concatenates gencode fasta with lrp fasta
# create a reference for each sample name from RNA analysis or custom count matrix
Rscript /scratch/cwp5au/LRP2_lite/bin/build_mass_spec_reference.R \
  --lrp_fasta /scratch/cwp5au/peptide_testing/data/test.predicted.proteome.all_best_orfs.fa \
  --gencode_version 47 \
  --sample_name NBM15 \
  --counts /scratch/cwp5au/peptide_testing/data/all_AML_all_hashids_with_cpm.txt \
  --outdir /scratch/cwp5au/peptide_testing/results/