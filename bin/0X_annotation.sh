#!/bin/bash

#SBATCH --job-name=02_annotation
#SBATCH --cpus-per-task=1          
#SBATCH --nodes=1                   
#SBATCH --ntasks-per-node=1         
#SBATCH --mem=32G           
#SBATCH --time=02:00:00             
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --output=log_files/%x-%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cwp5au@virginia.edu

set -e
source config/lrp2_config.sh
THREADS=8 # should match above
