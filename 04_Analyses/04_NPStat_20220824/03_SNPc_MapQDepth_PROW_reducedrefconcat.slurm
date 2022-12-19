#!/bin/bash

#SBATCH --job-name=03_SNPc_MapQDepth_PROW_reducedrefconcat
#SBATCH --partition=pcon06
#SBATCH --output=03_SNPc_MapQDepth_PROW_reducedrefconcat_%j.txt
#SBATCH --error=03_SNPc_MapQDepth_PROW_reducedrefconcat_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

# I created a samtools env with the newest version (1.15.1) using the following:

# ssh pinnacle-l6
# module load python/anaconda-3.9
# source /share/apps/bin/conda-3.9.sh
# mamba create -n samtools samtools

# so load the stuff needed to activate that environment below 

module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate samtools

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
# it's a version of this reference, best assembly from spades on 20211014 of PROW 981 R1: /scrfs/storage/amatthews/20210816_projects/20210816_exp/01_ASSEMBLY/spades_20211014/PROW_981_R1_TGCTGTGA-GAATCGTG/scaffolds.fasta.... but I have removed the microbe seqs and all contigs <5000bp (scaffolds_reduced_contigs_kept.fasta) and then concatenated all contigs into a super contig (scaffolds_reduced_contigs_kept_concatenated.fasta)

REF_DIR=$PROJECT_DIR/02_IndexRef/ref_full
REF=scaffolds_reduced_contigs_kept_concatenated
REFTYPE=ref_reduced_concat


## SNP DIR 
SNP_DIR=$PROJECT_DIR/03_SNP_20220513

## RESULTS DIR PER SPECIES
RESULTS_DIR=$SNP_DIR/RESULTS_${REFTYPE}_$SPP



##----------------------------------------------------------------------------------------
## Loop through samples
## Skip alignments with MAPQ smaller than 20, calculate average depth, and index *align_sort_dm_mq20.bams

while IFS=" " read -r value1 value2
do {

rsync -av $RESULTS_DIR/bam/${value1}_align_sort_dm.bam .

# skips alignments with mapping quality less than 20
samtools view --min-MQ 20 $RESULTS_DIR/bam/${value1}_align_sort_dm.bam --bam --output $RESULTS_DIR/bam/${value1}_align_sort_dm_mq20.bam

# calculate average depth and output to a .txt file
# -a means output all positions
samtools depth -a $RESULTS_DIR/bam/${value1}_align_sort_dm_mq20.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average =",sum/NR; print "Stdev =",sqrt(sumsq/NR - (sum/NR)**2)}' >> $RESULTS_DIR/bam/${value1}_align_sort_dm_mq20_depth.txt;

# index *_align_sort_dm_mq20.bam files (i.e., aligned, sorted, duplMarked, mapping quality of 20)
samtools index $RESULTS_DIR/bam/${value1}_align_sort_dm_mq20.bam



} done <"$filenames"






