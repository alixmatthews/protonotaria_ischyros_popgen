#!/bin/bash

#SBATCH --job-name=03_SNPe_mpileup_PROW_reducedref
#SBATCH --partition=comp72
#SBATCH --output=03_SNPe_mpileup_PROW_reducedref_%j.txt
#SBATCH --error=03_SNPe_mpileup_PROW_reducedref_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=10:00:00

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
## Do mpileup with supplying a list that has the absolute paths of the files
## -Ou flag indicates uncompressed bcf output to speed up performance by removing unnecessary compression/decompression 
## max-depth flag is related to: "at a position, read maximally INT reads per input file... Max raw per-file depth; avoids excessive memory usage". The default is 250, but they suggest the user has responsibility to change it based on their data. Our avg depth is around 15-16 for both spp. and suggest to 2-3x the average depth, so maybe 50 is good. You can also filter based on depth later (using GATK hard filters). Point is to avoid really high coverages (above the average) because that can suggest misalignment and higher changes of erroneous variant calls.

bcftools mpileup -Ou --fasta-ref ${REF_DIR}/${REF}.fasta --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --bam-list ${RESULTS_DIR}/bam/${SPP}_align_sort_dm_mq20_bam.list --max-depth 50 > ${RESULTS_DIR}/bam/${SPP}_ALL_align_sort_dm_mq20_bam.mpileup







