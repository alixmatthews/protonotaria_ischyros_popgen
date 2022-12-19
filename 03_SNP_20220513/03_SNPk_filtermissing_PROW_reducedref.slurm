#!/bin/bash

#SBATCH --job-name=03_SNPk_filtermissing_PROW_reducedref
#SBATCH --partition=comp01
#SBATCH --output=03_SNPk_filtermissing_PROW_reducedref_%j.txt
#SBATCH --error=03_SNPk_filtermissing_PROW_reducedref_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load vcftools/0.1.15

cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## "ADJUSTABLES"

## SPECIES OF INTEREST, PROW IN THIS CASE THE FILTERED DATASETS (REMOVED SAMPLES WITH HIGH PERCENTAGE OF POOR MAPPING)
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



##----------------------------------------------------------------------------------------
## Filter out SNP sites that have different levels of missing data (genotyping successful in 50%, 75%, and 100% of individuals... 50% would be least stringent (keeping more sites with missing data across individuals) versus 100% which would not allow any missing data for any individual)


# DP threshold of 5
vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp05.vcf.gz --max-missing 0.5 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss50

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp05.vcf.gz --max-missing 0.75 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss75

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp05.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100



# DP threshold of 10
vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10.vcf.gz --max-missing 0.5 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss50

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10.vcf.gz --max-missing 0.75 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss75

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100



# DP threshold of 15
vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15.vcf.gz --max-missing 0.5 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss50

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15.vcf.gz --max-missing 0.75 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss75

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100



# DP threshold of 20
vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp20.vcf.gz --max-missing 0.5 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss50

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp20.vcf.gz --max-missing 0.75 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss75

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp20.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100

