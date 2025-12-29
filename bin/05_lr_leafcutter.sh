#!/bin/bash

#SBATCH --job-name=05_lr_leafcutter
#SBATCH --cpus-per-task=1 #number of cores to use
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
#SBATCH --time=12:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --output=log_files/%x-%j.log

set -e
source config/lrp2_config.sh
THREADS=1 # should match above

module purge
module load gcc/11.4.0
module load openmpi/4.1.4
module load R/4.5.0
module load python/3.11.4

mkdir -p ${OUTPUT_DIR}/multisample_analysis/lr_leafcutter

#export R_LIBS_USER="" 
R -e ".libPaths()"
export R_LIBS_USER="/sfs/gpfs/tardis/home/cwp5au/R/goolf/4.5:$R_LIBS_USER"
R -e ".libPaths()"

# Step 1: Run long read leafcutter clustering (minicutter)
Rscript scripts/05_lr_leafcutter.R \
  --gtf ${OUTPUT_DIR}/sqanti_transcript/${OUTPUT_BASE_NAME}_corrected_filtered.gtf
  --counts ${OUTPUT_DIR}/sqanti_transcript/${OUTPUT_BASE_NAME}_hashids_with_cpm_filtered.txt
  --mode exon

# Step 2: Differential splicing
python -m venv environments/leafcutter_env
source environments/leafcutter_env/bin/activate
pip install -r environments/leafcutter_requirements.txt

python scripts/leafcutter_ds.py \
  -0 $CONTROL_GROUP \
  --min_samples_per_intron 2 \
  --min_samples_per_group 2 \
  --num_threads $THREADS \
  -o ${OUTPUT_DIR}/multisample_analysis/lr_leafcutter/${OUTPUT_BASE_NAME} \
  ${OUTPUT_DIR}/multisample_analysis/lr_leafcutter/${OUTPUT_BASE_NAME}_lr_leafcutter_perind_numers.counts.gz \
  ${OUTPUT_DIR}/multisample_analysis/lr_leafcutter/${OUTPUT_BASE_NAME}_groups_file.txt
  
deactivate