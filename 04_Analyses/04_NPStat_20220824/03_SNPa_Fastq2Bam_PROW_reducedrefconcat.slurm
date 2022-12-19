#!/bin/bash

#SBATCH --job-name=03_SNPa_Fastq2Bam_PROW_reducedrefconcat
#SBATCH --partition=pcon06
#SBATCH --output=03_SNPa_Fastq2Bam_PROW_reducedrefconcat_%j.txt
#SBATCH --error=03_SNPa_Fastq2Bam_PROW_reducedrefconcat_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load bwa/0.7.17
module load samtools/1.13
module load java/sunjdk_1.8.0 # for gatk and picard
module load gatk/4.2.6.1
module load fastqc/0.11.5 # workaround for picard-tools
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

## FASTQ LOCATIONS
FASTQ=$PROJECT_DIR/01_DataCombo_20210903/Adapter_Removed_bb

## REF DIR, REF, AND REFTYPE
# it's a version of this reference, best assembly from spades on 20211014 of PROW 981 R1: /scrfs/storage/amatthews/20210816_projects/20210816_exp/01_ASSEMBLY/spades_20211014/PROW_981_R1_TGCTGTGA-GAATCGTG/scaffolds.fasta.... but I have removed the microbe seqs and all contigs <5000bp (scaffolds_reduced_contigs_kept.fasta) and then concatenated all contigs into a super contig (scaffolds_reduced_contigs_kept_concatenated.fasta)

REF_DIR=$PROJECT_DIR/02_IndexRef/ref_full
REF=scaffolds_reduced_contigs_kept_concatenated
REFTYPE=ref_reduced_concat

## SNP DIR (mkdir'd beforehand)
SNP_DIR=$PROJECT_DIR/03_SNP_20220513

## RESULTS DIR PER SPECIES
RESULTS_DIR=$SNP_DIR/RESULTS_${REFTYPE}_$SPP

## MAKE SUBDIRS IN RESULTS_DIR
mkdir -p $RESULTS_DIR # necessary to make multiple dirs within
mkdir $RESULTS_DIR/bam
mkdir $RESULTS_DIR/sam
mkdir $RESULTS_DIR/flagstat




##----------------------------------------------------------------------------------------
## LOOP THROUGH SAMPLES
## fastq to indexed, align_sort.bam

while IFS=" " read -r value1 value2
do {

# align fastq reads to reference genome using BWA, output is a .sam file
# -M: tells bwa to consider split reads as secondary, required for GATK variant calling
# -R: read group info
bwa mem -t 32 -M -R "@RG\tID:${value1}\tSM:${value1}\tLB:${value2}\tPL:ILLUMINA" $REF_DIR/$REF.fasta $FASTQ/${value1}/${value1}_bb_trim_R1.fastq $FASTQ/${value1}/${value1}_bb_trim_R2.fastq > ${value1}_align.sam

# convert .sam to .bam file (bam is a compressed binary version of sam) using samtools
samtools view --threads 30 -S -b ${value1}_align.sam > ${value1}_align.bam

# sort the bam files using samtools
samtools sort --threads 30 ${value1}_align.bam -o ${value1}_align_sort.bam 

# index the align_sort.bam files
samtools index ${value1}_align_sort.bam

# mapping stats
samtools flagstat --threads 30 -O tsv ${value1}_align_sort.bam > ${value1}_align_sort_flagstat_out.txt


# rsync results back over to home directory
rsync -av --remove-source-files ${value1}_align.sam $RESULTS_DIR/sam/
rsync -av --remove-source-files ${value1}_align.bam $RESULTS_DIR/bam/
rsync -av --remove-source-files ${value1}_align_sort.bam $RESULTS_DIR/bam/
rsync -av --remove-source-files ${value1}_align_sort_flagstat_out.txt $RESULTS_DIR/flagstat/


} done <"$filenames"


