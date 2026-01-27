#!/bin/bash

#SBATCH --job-name=21_novel_peptides
#SBATCH --cpus-per-task=10 #number of cores to use
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --account=sheynkman_lab
#SBATCH --output=%x-%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yqy3cu@virginia.edu

module load gcc/11.4.0  
module load openmpi/4.1.4
module load python/3.11.4
module load bioconda/py3.10
module load miniforge/24.3.0-py3.11

conda activate reference_tab

# Refined
python ./00_scripts/17_novel_peptides.py \
--pacbio_peptides ./16_MetaMorpheus/AllPeptides.jurkat.refined.psmtsv \
--gencode_fasta ./02_make_gencode_database/gencode_clusters.fasta \
--uniprot_fasta ./00_input_data/uniprot_reviewed_canonical_and_isoform.fasta \
--name ./17_novel_peptides/jurkat_refined

# Filtered
python ./00_scripts/17_novel_peptides.py \
--pacbio_peptides ./16_MetaMorpheus/AllPeptides.jurkat.filtered.psmtsv \
--gencode_fasta ./02_make_gencode_database/gencode_clusters.fasta \
--uniprot_fasta ./00_input_data/uniprot_reviewed_canonical_and_isoform.fasta \
--name ./17_novel_peptides/jurkat_filtered

# Hybrid
python ./00_scripts/21_novel_peptides.py \
--pacbio_peptides ./16_MetaMorpheus/AllPeptides.jurkat.hybrid.psmtsv \
--gencode_fasta ./02_make_gencode_database/gencode_clusters.fasta \
--uniprot_fasta ./00_input_data/uniprot_reviewed_canonical_and_isoform.fasta \
--name ./17_novel_peptides/jurkat_hybrid

conda deactivate
module purge
