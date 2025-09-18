#!/bin/bash

#SBATCH --job-name=03_orf_calling
#SBATCH --cpus-per-task=1 #number of cores to use
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00 #amount of time for the whole job
#SBATCH --mem=16G
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --output=log_files/%x-%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cwp5au@virginia.edu #your email address to receive notifications

set -e
source config/lrp2_config.sh
THREADS=1 # should match above

module purge
module load gcc/11.4.0
module load openmpi/4.1.4
module load miniforge/24.3.0-py3.11

mkdir -p ${OUTPUT_DIR}/orf_calling
source $(conda info --base)/etc/profile.d/conda.sh

# Check if environment exists, create if needed. This can take awhile.
if ! conda info --envs | grep -q "$LRP2_ENV_NAME"; then
    echo "Creating conda environment $LRP2_ENV_NAME from $LRP2_ENV_FILE..."
    conda env create -f "$LRP2_ENV_FILE" -n "$LRP2_ENV_NAME"
else
    echo "Environment $LRP2_ENV_NAME already exists."
fi

conda activate "$LRP2_ENV_NAME"
export R_LIBS_USER="" 
R -e ".libPaths()"

# Use specific files based on species
case $SPECIES in
    "human")
        echo "Using Human Hexamer file and logit model..."
        HEXAMER_INPUT="external_data/Human_Hexamer.tsv"
        MODEL_INPUT="external_data/Human_logitModel.RData"
        export CPAT_CODING_THRESHOLD=0.364
        ;;
    "mouse")
        echo "Using Mouse Hexamer file and logit model..."
        HEXAMER_INPUT="external_data/Mouse_Hexamer.tsv"
        MODEL_INPUT="external_data/Mouse_logitModel.RData"
        export CPAT_CODING_THRESHOLD=0.44
        ;;
esac

# Step 1: Run CPAT
cpat.py \
   -x $HEXAMER_INPUT \
   -d $MODEL_INPUT \
   -g ${OUTPUT_DIR}/sqanti/${OUTPUT_BASE_NAME}_corrected_filtered.fasta \
   --min-orf=$MIN_ORF \
   --top-orf=$TOP_ORF \
   -o ${OUTPUT_DIR}/orf_calling/${OUTPUT_BASE_NAME}_cpat \
   2> ${OUTPUT_DIR}/orf_calling/${OUTPUT_BASE_NAME}_cpat.error 

# Step 2: Filter CPAT output to call best ORF 
Rscript 03_filter_cpat.R

conda deactivate