#!/bin/bash

#SBATCH --job-name=00_PP_SNPDATA_multiqc_before_q_trim
#SBATCH --partition=pcon06
#SBATCH --output=00_PP_SNPDATA_multiqc_before_q_trim_%j.txt
#SBATCH --error=00_PP_SNPDATA_multiqc_before_q_trim_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00

module load java/sunjdk_1.8.0
module load fastqc/0.11.5

module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate bioconda3-el7

cd /local_scratch/$SLURM_JOB_ID/


############ MULTIQC ############

## run multiqc - just give it the path to fastqc files
multiqc /scrfs/storage/amatthews/20210816_projects/20210816_snp/00_PP_20210902/FASTQC/Before_Q_Trim

multiqc /scrfs/storage/amatthews/20210412_snp/00_PP_20210427/FASTQC/Before_Q_Trim


