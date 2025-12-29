#!/bin/bash

#SBATCH --job-name=16_metamorph
#SBATCH --cpus-per-task=10 #number of cores to use
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00 #amount of time for the whole job
#SBATCH --partition=standard #the queue/partition to run on
#SBATCH --account=sheynkman_lab
#SBATCH --output=%x-%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yqy3cu@virginia.edu

module load gcc/11.4.0  
module load openmpi/4.1.4
module load python/3.11.4
module load bioconda/py3.10
module load anaconda/2023.07-py3.11

conda activate metamorph

#gencode only
metamorpheus -t ./00_input_data/Task1SearchTaskconfig_orf.toml \
-s ./15_MS_file_convert/120426_Jurkat_highLC_Frac1.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac2.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac3.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac4.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac5.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac6.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac7.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac8.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac9.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac10.mzML \
-d ./02_make_gencode_database/gencode_clusters.fasta \
-o ./16_MetaMorpheus/gencode/

#hybrid database
metamorpheus -t ./00_input_data/Task1SearchTaskconfig_orf.toml \
-s ./15_MS_file_convert/120426_Jurkat_highLC_Frac1.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac2.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac3.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac4.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac5.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac6.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac7.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac8.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac9.mzML ./15_MS_file_convert/120426_Jurkat_highLC_Frac10.mzML \
-d ./13_protein_hybrid_database/jurkat_hybrid.fasta \
-o ./16_MetaMorpheus/hybrid/

conda deactivate
