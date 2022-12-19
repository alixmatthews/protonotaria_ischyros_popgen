#!/bin/bash

#SBATCH --job-name=04_NPStat_pileupbysample_PROW_reducedrefconcat
#SBATCH --partition=pcon06
#SBATCH --output=04_NPStat_pileupbysample_PROW_reducedrefconcat_%j.txt
#SBATCH --error=04_NPStat_pileupbysample_PROW_reducedrefconcat_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load samtools/1.15.1

cd /local_scratch/$SLURM_JOB_ID/

##----------------------------------------------------------------------------------------
## "ADJUSTABLES"

## SPECIES OF INTEREST
SPP=PROW 

## SAMPLE NAMES - ABSOLUTE PATH
samples="/scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/PROW_FILTERED_samples_hapsize_NPStat.txt"

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_snp

## REF DIR, REF, AND REFTYPE
# it's a version of this reference, best assembly from spades on 20211014 of PROW 981 R1: /scrfs/storage/amatthews/20210816_projects/20210816_exp/01_ASSEMBLY/spades_20211014/PROW_981_R1_TGCTGTGA-GAATCGTG/scaffolds.fasta.... but I have removed the microbe seqs and all contigs <5000bp (scaffolds_reduced_contigs_kept.fasta) and then concatenated all contigs into a super contig (scaffolds_reduced_contigs_kept_concatenated.fasta)

REF_DIR=$PROJECT_DIR/02_IndexRef/ref_full
REF=scaffolds_reduced_contigs_kept_concatenated
REFTYPE=ref_reduced_concat

## SNP DIR 
SNP_DIR=$PROJECT_DIR/03_SNP_20220513

## RESULTS DIR PER SPECIES
RESULTS_DIR=$SNP_DIR/RESULTS_${REFTYPE}_$SPP



##----------------------------------------------------------------------------------------
## make pileup files of each sample individually


cd $RESULTS_DIR/bam

while IFS=" " read -r sample_name haploidy;
do samtools mpileup ${sample_name}_align_sort_dm_mq20.bam > ${sample_name}_align_sort_dm_mq20.pileup; 
done<"$samples"




