#!/bin/bash

#SBATCH --job-name=04_plink_LDvari_PCA_PROW_reducedref
#SBATCH --partition=comp01
#SBATCH --output=04_plink_LDvari_PCA_PROW_reducedref_%j.txt
#SBATCH --error=04_plink_LDvari_PCA_PROW_reducedref_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load plink

# version is PLINK v1.90b5.2

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
mkdir ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20220706
PLINK_DIR=${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20220706



##----------------------------------------------------------------------------------------
## Use PLINK to identify variants that are potentially in linkage disequilibrium and separate them. Output files will be: the prune.in (variants to keep) and prune.out (variants that are to be removed because potentially linked)

## Also make a .bim and .bed file for ADMIXTURE and .eigenval and .eigenvec files for plink-based PCA


## First identify variants that are potentially linked. Consider a window size of 100 SNPs, calculate LD between each pair of SNPs in the window, remove one of a pair of SNPs if the r2 value is greater than 0.001, 0.01, 0.1, 0.3, 0.4, then shift the window down 50 SNPs forward and repeat the procedure. 

## Do for 10X and 15X datasets (0% missing) so can compare PCAs


cd ${PLINK_DIR}


plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.001 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld0001

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.01 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.1 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld01

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.3 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld03

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.4 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld04



plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.001 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld0001

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.01 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.1 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.3 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld03

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.4 --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld04






## Then, make the .bim and .bed files for ADMIXTURE and files to create a PLINK-based PCA


plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld0001.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld0001

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld01.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld01

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld03.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld03

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld04.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld04



plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld0001.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld0001

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld03.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld03

plink --vcf ../${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld04.prune.in --make-bed --pca --out ${SPP}_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld04
