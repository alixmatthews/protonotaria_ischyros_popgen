#!/bin/bash

#SBATCH --job-name=02_IndexRef_concat_quast
#SBATCH --partition=pcon06
#SBATCH --output=02_IndexRef_concat_quast_%j.txt
#SBATCH --error=02_IndexRef_concat_quast_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:15:00


module load java/sunjdk_1.8.0
module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate bioconda3-el7
module load quast/5.0.0 # I think this is the only one I really need, but loading the others first just in case


cd /local_scratch/$SLURM_JOB_ID/

### MAKE OUTPUT DIR ###
OUTPUT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_snp/02_IndexRef/ref_full/quast_compare_concatenated_contigs
mkdir ${OUTPUT_DIR}

### RUN QUAST ###


/share/apps/bioinformatics/quast/quast-5.0.0/quast.py --eukaryote -t 4 -o ${OUTPUT_DIR} /scrfs/storage/amatthews/20210816_projects/20210816_snp/02_IndexRef/ref_full/scaffolds_reduced_contigs_kept.fasta /scrfs/storage/amatthews/20210816_projects/20210816_snp/02_IndexRef/ref_full/scaffolds_reduced_contigs_kept_concatenated.fasta 

