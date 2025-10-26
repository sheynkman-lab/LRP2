#!/bin/bash

#SBATCH --job-name=01_isoseq
#SBATCH --cpus-per-task=8          
#SBATCH --nodes=1                   
#SBATCH --ntasks-per-node=1         
#SBATCH --mem=32G           
#SBATCH --time=12:00:00             
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --output=log_files/%x-%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cwp5au@virginia.edu

set -e
source config/lrp2_jurkat_test_config.sh
THREADS=8 # should match above

module purge
module load isoseqenv # pre-loaded isoseq for clustering
module load smrtlink

mkdir -p ${OUTPUT_DIR}/pacbio_isoseq # create one output directory, if not already created

# input directory and number of samples are specified in the config file 
INPUT_FILES=($(tail -n +2 "$SAMPLE_METADATA" | awk '{print $1}'))

echo "Starting Iso-Seq pipeline at $(date)"

# Step 1: Merge
echo "Merging ${#INPUT_FILES[@]} BAM files at $(date)..."
pbmerge -o ${OUTPUT_DIR}/pacbio_isoseq/${OUTPUT_BASE_NAME}_merged.flnc.bam "${INPUT_FILES[@]}"

# Step 2: Cluster (needs a lot of memory)
echo "Clustering reads at $(date)..."
isoseq cluster2 \
  ${OUTPUT_DIR}/pacbio_isoseq/${OUTPUT_BASE_NAME}_merged.flnc.bam \
  ${OUTPUT_DIR}/pacbio_isoseq/${OUTPUT_BASE_NAME}_merged.clustered.bam

# Step 3: Align
# Figure out how to get --split-by-sample, Biosample naming is maintained in merged flnc but not merged clustered bam
echo "Aligning to genome at $(date)..."
pbmm2 align --preset ISOSEQ --sort -j $THREADS --log-level INFO \
  $GENCODE_GENOME_FA \
  ${OUTPUT_DIR}/pacbio_isoseq/${OUTPUT_BASE_NAME}_merged.clustered.bam \
  ${OUTPUT_DIR}/pacbio_isoseq/${OUTPUT_BASE_NAME}_merged.aligned.bam

# Step 4: Collapse
echo "Collapsing reads at $(date)..."
isoseq collapse --max-fuzzy-junction $MAX_FUZZY_JUNCTION --max-5p-diff $MAX_5P_DIFF --max-3p-diff $MAX_3P_DIFF \
    ${OUTPUT_DIR}/pacbio_isoseq/${OUTPUT_BASE_NAME}_merged.aligned.bam \
    ${OUTPUT_DIR}/pacbio_isoseq/${OUTPUT_BASE_NAME}_merged.flnc.bam \
    ${OUTPUT_DIR}/pacbio_isoseq/${OUTPUT_BASE_NAME}_merged.collapsed.gff

echo "Iso-Seq pipeline completed at $(date)"
echo "Outputs saved to: $OUTPUT_DIR/pacbio_isoseq"
