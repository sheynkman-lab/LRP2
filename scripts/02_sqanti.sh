#!/bin/bash

#SBATCH --job-name=02_sqanti
#SBATCH --cpus-per-task=8 #number of cores to use
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00 #amount of time for the whole job
#SBATCH --mem=32G
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --output=log_files/%x-%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cwp5au@virginia.edu #your email address to receive notifications

set -e
source config/lrp2_config.sh
THREADS=8 # should match above

module purge
module load gcc/11.4.0
module load openmpi/4.1.4
module load miniforge/24.3.0-py3.11

mkdir -p ${OUTPUT_DIR}/sqanti
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

# get gtf based on source
case $GTF_SOURCE in
    "lrp_isoseq")
        echo "Using LRP Isoseq generated GTF..."
        LRS_GTF="results/pacbio_isoseq/merged.collapsed.gff"
        ;;
    "custom")
        echo "Using custom GTF: $CUSTOM_GTF_PATH"
        LRS_GTF="$CUSTOM_GTF_PATH"
        ;;
esac

# Step 1: Run SQANTI3 on long read gtf
python $SQANTI_PATH/sqanti3_qc.py \
    --force_id_ignore \
    --skipORF \
    --output $OUTPUT_BASE_NAME \
    --dir ${OUTPUT_DIR}/sqanti \
    --cpus $THREADS \
    --report both \
    --fl_count results/pacbio_isoseq/merged.collapsed.flnc_count.txt \
    --isoforms $LRS_GTF \
    --refGTF $GENCODE_GTF_FILE \
    --refFasta $GENCODE_GENOME_FA

# Step 2: Filter SQANTI transcripts with custom script
Rscript scripts/02_filter_sqanti_transcripts.R 

conda deactivate