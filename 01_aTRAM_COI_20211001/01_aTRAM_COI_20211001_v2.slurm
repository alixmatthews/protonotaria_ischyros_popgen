#!/bin/bash

#SBATCH --job-name=01_aTRAM_COI_20211001_v2
#SBATCH --partition=comp72
#SBATCH --output=01_aTRAM_COI_20211001_v2_%j.txt
#SBATCH --error=01_aTRAM_COI_20211001_v2_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00

module load exonerate/2.4.0
module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate aTRAM

cd /local_scratch/$SLURM_JOB_ID/

############################################################
############ CONFIGURATION INFORMATION ############
############ LOCATIONS AND SAMPLE NAMES ############
############################################################

## SAMPLE FILE NAMES 
# this file has all the CERW and PROW samples from SNP analyses, no second column with RG info. This file is every other sample from the previous file (all_cerw_and_prow_filenames.txt)... something weird happened with exonerate so I am re-running it with only the samples that did not get exonerated. Commenting out everything else to hopefully streamline the pipeline.

filenames="/scrfs/storage/amatthews/20210816_projects/20210816_phylo/all_cerw_and_prow_filenames_v2.txt" 

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_phylo

## ATRAM DIR
## contains COI reference (coi_Proctophyllodes_miliariae)
ATRAM_DIR=$PROJECT_DIR/01_aTRAM_COI_20211001

## ATRAM DATABASES
ATRAM_DB=$ATRAM_DIR/aTRAM_db
#mkdir ATRAM_DB


############ SET UP PER SAMPLE DIRS ############

while IFS=" " read -r value1
do {


## FASTQ FILE DIRECTORY
SAMPLE_FASTQ=/scrfs/storage/amatthews/20210816_projects/20210816_snp/01_DataCombo_20210903/Adapter_Removed_bb/${value1}


# ATRAM_DB DIR
#mkdir $ATRAM_DB/{value1}

# SAMPLE OUTPUT / COI / (ALL, FILTERED, STITCH)
SAMPLE_OUT=$ATRAM_DIR/${value1}_out
#mkdir $SAMPLE_OUT
SAMPLE_OUT_COI=$SAMPLE_OUT/${value1}_coi
#mkdir -p $SAMPLE_OUT_COI
#mkdir $SAMPLE_OUT_COI/all
#mkdir $SAMPLE_OUT_COI/filtered
#mkdir $SAMPLE_OUT_COI/stitch



############################################################
############ RUN aTRAM ############
############################################################

############ PREPROCESSOR ############

#atram_preprocessor.py --blast-db=$ATRAM_DB/${value1}/${value1} --end-1=$SAMPLE_FASTQ/${value1}_bb_trim_R1.fastq  --end-2=$SAMPLE_FASTQ/${value1}_bb_trim_R2.fastq
## second /${value1} on --blast-db line is referring to the prefix to use for the file names it creates in this step


############ ASSEMBLER ############
############ COI ############

#atram.py --blast-db=$ATRAM_DB/${value1}/${value1} --query-split=$ATRAM_DIR/COI_ref.fasta --assembler=spades --protein --db-gencode=5 --output-prefix=$SAMPLE_OUT_COI/all/
## second /${value1} on --blast-db line is referring to the prefix to use/search for the file names it created in the previous step

cd $SAMPLE_OUT_COI
#find . -name "*filtered_contigs.fasta" -exec mv -i {} -t "./filtered/" ";"
## move filtered_contigs files to a new folder called "filtered" (separate them from non-filtered)



############ STITCHER ############
############ COI ############

atram_stitcher.py --assemblies-dir=$SAMPLE_OUT_COI/filtered --reference-genes=$ATRAM_DIR/COI_ref.fasta --taxa=/scrfs/storage/amatthews/20210816_projects/20210816_phylo/all_cerw_and_prow_filenames_v2.txt --output-prefix=$SAMPLE_OUT_COI/stitch/${value1}



############ FINISH aTRAM ############

} done <"$filenames"

