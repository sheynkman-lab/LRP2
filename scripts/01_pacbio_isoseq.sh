#!/bin/bash

#SBATCH --job-name=01_isoseq
#SBATCH --cpus-per-task=40          
#SBATCH --nodes=1                   
#SBATCH --ntasks-per-node=1         
#SBATCH --mem=1460G           
#SBATCH --time=72:00:00             
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --output=log_files/%x-%j.log

set -e
source config/lrp2_config.sh
THREADS=40 # should match above

module purge
module load isoseqenv # pre-loaded isoseq for clustering
module load smrtlink
module load samtools

mkdir -p ${OUTPUT_DIR}/pacbio_isoseq # create one output directory, if not already created

# input directory and number of samples are specified in the config file 
INPUT_FILES=($(tail -n +2 "$SAMPLE_METADATA" | awk '{print $1}'))

echo "Starting Iso-Seq pipeline at $(date)"

# Step 1: Merge
echo "Merging ${#INPUT_FILES[@]} BAM files..."
pbmerge -o ${OUTPUT_DIR}/pacbio_isoseq/merged.flnc.bam "${INPUT_FILES[@]}"

# Step 2: Cluster
echo "Clustering reads..."
isoseq cluster2 ${OUTPUT_DIR}/pacbio_isoseq/merged.flnc.bam ${OUTPUT_DIR}/pacbio_isoseq/merged.clustered.bam

# Step 3: Align
echo "Aligning to genome..."
pbmm2 align $GENOME_FA ${OUTPUT_DIR}/pacbio_isoseq/merged.clustered.bam ${OUTPUT_DIR}/pacbio_isoseq/merged.aligned.bam --preset ISOSEQ --split-by-sample --sort -j $THREADS --log-level INFO

# Step 4: Collapse
echo "Collapsing reads..."
isoseq collapse --max-fuzzy-junction $MAX_FUZZY_JUNCTION --max-5p-diff $MAX_5P_DIFF --max-3p-diff $MAX_3P_DIFF \
    ${OUTPUT_DIR}/pacbio_isoseq/merged.aligned.bam \
    ${OUTPUT_DIR}/pacbio_isoseq/merged.flnc.bam \
    ${OUTPUT_DIR}/pacbio_isoseq/merged.collapsed.gff

echo "Iso-Seq pipeline completed at $(date)"
echo "Outputs saved to: $OUTPUT_DIR/pacbio_isoseq"
