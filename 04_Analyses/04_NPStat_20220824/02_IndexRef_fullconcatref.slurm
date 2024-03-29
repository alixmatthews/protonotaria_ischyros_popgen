#!/bin/bash

#SBATCH --job-name=02_IndexRef_fullconcatref
#SBATCH --partition=pcon06
#SBATCH --output=02_IndexRef_fullconcatref_%j.txt
#SBATCH --error=02_IndexRef_fullconcatref_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

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
## MAKE AND ASSIGN DIRECTORIES

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_snp

## SNP DIR 
SNP_DIR=$PROJECT_DIR/02_IndexRef

## REF DIR (new subdir for the full genome)
REF_DIR=$SNP_DIR/ref_full

## REFERENCE GENOME
REF_FULL=scaffolds_reduced_contigs_kept_concatenated


##----------------------------------------------------------------------------------------
## INDEX, .FAI, AND .DICT THE REFERENCES
## only have to do this once per reference
## REF_FULL
bwa index $REF_DIR/$REF_FULL.fasta

## ensure .fai file is created
samtools faidx $REF_DIR/$REF_FULL.fasta

## create sequence dictionary for reference
java -jar /share/apps/bioinformatics/picard/picard-tools-2.17.10/picard.jar CreateSequenceDictionary REFERENCE=$REF_DIR/$REF_FULL.fasta OUTPUT=$REF_DIR/$REF_FULL.dict



