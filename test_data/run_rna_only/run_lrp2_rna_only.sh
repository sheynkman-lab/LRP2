#!/bin/bash

#SBATCH --job-name=nf_lrp2_driver
#SBATCH --cpus-per-task=1
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH --mem=20GB
#SBATCH --output=log_files/slurm-%j.out

module load apptainer nextflow

nextflow run main.nf \
    --input test_data/run_rna_only/samplesheet_rna_only.csv \
    --outdir test_data/run_rna_only/test_results \
    --dataset_name lrptest \
    --genome GRCh38.p14.v49 \
    -profile singularity,slurm
    