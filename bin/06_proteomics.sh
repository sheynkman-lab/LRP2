#!/bin/bash

#SBATCH --job-name=06_proteomics
#SBATCH --cpus-per-task=10 #number of cores to use
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --account=sheynkman_lab
#SBATCH --output=%x-%j.log

module load gcc/11.4.0  
module load openmpi/4.1.4
module load python/3.11.4
module load bioconda/py3.10
module load miniforge/24.3.0-py3.11

set -e
source config/lrp2_config.sh
THREADS=1 # should match above

mkdir -p ${OUTPUT_DIR}/proteomics

conda activate metamorph
conda install -c conda-forge metamorpheus

metamorpheus -t config/Task1SearchTaskconfig_orf.toml \
-s data/120426_Jurkat_highLC_Frac1.mzML data/120426_Jurkat_highLC_Frac2.mzML \
-d ${OUTPUT_DIR}/proteomics/hybrid_database.fasta \
-o ${OUTPUT_DIR}/proteomics/

conda deactivate