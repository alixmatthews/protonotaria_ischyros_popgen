#!/bin/bash

#SBATCH --job-name=03_SNPb_MrkDup_PROW_reducedref
#SBATCH --partition=comp06
#SBATCH --output=03_SNPb_MrkDup_PROW_reducedref_%j.txt
#SBATCH --error=03_SNPb_MrkDup_PROW_reducedref_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load samtools/1.13
module load java/sunjdk_1.8.0 # for gatk and picard
module load gatk/4.2.6.1
module load fastqc/0.11.5 # workaround for picard-tools to load
module load picard-tools/2.17.10

cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## "ADJUSTABLES"

## SPECIES OF INTEREST
SPP=PROW 

## FILENAMES
filenames="/scrfs/storage/amatthews/20210816_projects/20210816_snp/all_prow_filenames.txt" 

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
## Loop through samples
## mark duplicates, index *align_sort_dm.bams

while IFS=" " read -r value1 value2
do {

rsync -av $RESULTS_DIR/bam/${value1}_align_sort.bam .

# mark duplicates in aligned_sorted.bam files
# this marks duplicates so that variant calling tools will ignore them. They are not removed from the dataset unless requested.
java -jar /share/apps/bioinformatics/picard/picard-tools-2.17.10/picard.jar MarkDuplicates INPUT=$RESULTS_DIR/bam/${value1}_align_sort.bam OUTPUT=$RESULTS_DIR/bam/${value1}_align_sort_dm.bam METRICS_FILE=$RESULTS_DIR/bam/${value1}_align_sort_dm_metrics.txt

# index *_align_sort_dm.bam files (i.e., aligned, sorted, duplMarked)
samtools index $RESULTS_DIR/bam/${value1}_align_sort_dm.bam


} done <"$filenames"






