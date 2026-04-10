#!/bin/bash

#SBATCH --job-name=nf_lrp2_driver
#SBATCH --cpus-per-task=1
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH --mem=20GB
#SBATCH --output=log_files/slurm-%j.out

module load apptainer nextflow

nextflow run main.nf \
    --input test_data/run_fragpipe_rna/samplesheet_fragpipe_rna.csv \
    --outdir test_data/run_fragpipe_rna/test_results \
    --dataset_name lrptest \
    --genome GRCh38.p14.v49 \
    --protein_search fragpipe \
    --fragpipe_first_name "Megan" \
    --fragpipe_last_name "Schertzer" \
    --fragpipe_email "cwp5au@virginia.edu" \
    --fragpipe_institution "University of Virginia" \
    --fragpipe_token "101690" \
    --fragpipe_license_accept true \
    -profile singularity,slurm
    