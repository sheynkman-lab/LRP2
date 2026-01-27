#!/bin/bash

#SBATCH --job-name=env_test
#SBATCH --time=04:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --mem=32G
#SBATCH --output=log_files/%x-%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cwp5au@virginia.edu #your email address to receive notifications

module purge
module load gcc/11.4.0
module load openmpi/4.1.4

export PATH="$HOME/miniforge3/bin:$PATH"
source "$HOME/miniforge3/etc/profile.d/conda.sh"

which mamba
mamba --version

mamba env create -f environments/LRP2_update.conda_env.yml -n LRP2_update