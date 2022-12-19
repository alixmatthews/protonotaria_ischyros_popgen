#!/bin/bash

#SBATCH --job-name=04_ADMIXTURE_PROW_15x100missing_ld01_reducedref
#SBATCH --partition=comp01
#SBATCH --output=04_ADMIXTURE_PROW_15x100missing_ld01_reducedref_%j.txt
#SBATCH --error=04_ADMIXTURE_PROW_15x100missing_ld01_reducedref_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

## DON'T BELIEVE I NEED ANY FOR THIS SCRIPT

cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## "ADJUSTABLES"

## SPECIES OF INTEREST
SPP=PROW 

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


## PLINK DIR PER SPECIES
PLINK_DIR=${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20220706


## ADMIXTURE DIR PER SPECIES
mkdir ${RESULTS_DIR}/vcf/${SPP}_ADMIXTURE_20220709
ADMIXTURE_DIR=${RESULTS_DIR}/vcf/${SPP}_ADMIXTURE_20220709



##----------------------------------------------------------------------------------------
## Run ADMIXTURE with 10x cross-validation for K = 1-10 on 15X datasets with LD01


## First, do some reformatting of the .bim and rename as the original (to use for ADMIXTURE). This is because 'ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0 (found from https://speciationgenomics.github.io/ADMIXTURE/)'... Go ahead and do for all datasets at this point.


cd ${PLINK_DIR}

awk '{$1=0;print $0}' ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.bim > ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.bim.tmp
mv ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.bim.tmp ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.bim






## Now, run ADMIXTURE with 10x cross-validation for K = 1-10 on 15X datasets with LD01


cd ${ADMIXTURE_DIR}


for K in {1..10}
do
 /share/apps/bioinformatics/admixture/1.3.0/admixture --cv=10 ${PLINK_DIR}/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.bed $K > log_${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01_${K}.out
done







