#!/bin/bash

#SBATCH --job-name=05_multisample_analysis
#SBATCH --cpus-per-task=1 #number of cores to use
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
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

mkdir -p ${OUTPUT_DIR}/multisample_analysis

#export R_LIBS_USER="" 
R -e ".libPaths()"
export R_LIBS_USER="/sfs/gpfs/tardis/home/cwp5au/R/goolf/4.5:$R_LIBS_USER"
R -e ".libPaths()"

# Step 1: Run edgeR and DRIMSeq differential analysis
Rscript scripts/05_multisample_analysis.R