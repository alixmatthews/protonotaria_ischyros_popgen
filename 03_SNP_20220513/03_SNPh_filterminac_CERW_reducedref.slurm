#!/bin/bash

#SBATCH --job-name=03_SNPh_filterminac_CERW_reducedref
#SBATCH --partition=comp01
#SBATCH --output=03_SNPh_filterminac_CERW_reducedref_%j.txt
#SBATCH --error=03_SNPh_filterminac_CERW_reducedref_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

# I created a bcftools env with the newest version (1.15.1) using the following:

# ssh pinnacle-l6
# module load python/anaconda-3.9
# source /share/apps/bin/conda-3.9.sh
# mamba create -n BCFTools bcftools

# so load the stuff needed to activate that environment below 

module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate BCFTools

cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## "ADJUSTABLES"

## SPECIES OF INTEREST
SPP=CERW 

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_snp

## REF DIR, REF, AND REFTYPE
# it's a version of this reference, best assembly from spades on 20211014 of PROW 981 R1: /scrfs/storage/amatthews/20210816_projects/20210816_exp/01_ASSEMBLY/spades_20211014/PROW_981_R1_TGCTGTGA-GAATCGTG/scaffolds.fasta.... but I have removed the microbe seqs and all contigs <5000bp (scaffolds_reduced_contigs_kept.fasta)

REF_DIR=$PROJECT_DIR/02_IndexRef/ref_full
REF=scaffolds_reduced_contigs_kept
REFTYPE=ref_reduced


## SNP DIR 
SNP_DIR=$PROJECT_DIR/03_SNP_20220513

## RESULTS DIR PER SPECIES
RESULTS_DIR=$SNP_DIR/RESULTS_${REFTYPE}_$SPP



##----------------------------------------------------------------------------------------
## Filter out SNPs that different from the reference but are equal to each other across samples

bcftools view --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1.vcf.gz --min-ac 1:minor ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30.vcf.gz








