#!/bin/bash

#SBATCH --job-name=04_sqanti_protein
#SBATCH --cpus-per-task=8 #number of cores to use
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --output=log_files/%x-%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cwp5au@virginia.edu

set -e
source config/lrp2_config.sh
THREADS=8 # should match above

module purge
module load gcc/11.4.0
module load openmpi/4.1.4
module load miniforge/24.3.0-py3.11

mkdir -p ${OUTPUT_DIR}/protein_sqanti

# Check if environment exists, create if needed. This can take awhile.
if ! conda info --envs | grep -q "$LRP2_ENV_NAME"; then
    echo "Creating conda environment $LRP2_ENV_NAME from $LRP2_ENV_FILE..."
    mamba env create -f "$LRP2_ENV_FILE" -n "$LRP2_ENV_NAME"
else
    echo "Environment $LRP2_ENV_NAME already exists."
fi

source $(conda info --base)/etc/profile.d/conda.sh
mamba activate "$LRP2_ENV_NAME"
unset R_LIBS_USER

#export R_LIBS_USER="" 
R -e ".libPaths()"

# python scripts/sqanti3_protein_input_full_gtf.py \
#   ${OUTPUT_DIR}/orf_calling/${OUTPUT_BASE_NAME}_corrected_filtered_CDS.gtf \
#   ${OUTPUT_DIR}/orf_calling/${OUTPUT_BASE_NAME}_best_orfs_mapped.tsv \
#   $GENCODE_GTF_FILE \
#   -d ${OUTPUT_DIR}/protein_sqanti \
#   -p $OUTPUT_BASE_NAME

mamba deactivate