#!/bin/bash

#SBATCH --job-name=novel_peptides
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

Rscript /scratch/cwp5au/LRP2_lite/bin/novel_peptides.R \
  --sample_name A549 \
  --ms_search_software fragpipe \
  --acquisition_type DIA \
  --outdir /scratch/cwp5au/ms_test