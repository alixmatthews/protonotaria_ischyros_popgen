## NPStat analysis pipeline

- NPStat requires files to be a single population (sample) and from a single chromosome/contig/scaffold. 

- HOWEVER, the 20220809 version of my NPStat analysis is not quite what we need to do. See notes below.

- Luca Ferretti (NPStat developer) said that NPStat is not really well suited for contigs. His suggestion was to create a linear "pseudochromosome" of all the contigs concatenated to one another. Then align the reads to this pseudochromosome, treat it as a single chromosome, then feed the mpileup file(s) of each sample to NPStat. He also said that if we have a filtered VCF and you would like to use only those sites (which is true in our case), do not extract them from the pileup because NPStat needs the full pileup. Instead, we need to save the positions of all the sites into a file (one per each row) and feed it to NPStat using the "-snpfile" option. We have to convert the .bed file to consider all sites as unique positions and that there is only a single chromosome to make this -snpfile (right now site positions are not necessarily unique/may be repeated because all sites are separated by contig in the .bed file). 

- So let's get started!

--- 
### 1. Create a linear pseudochromosome

- *Concept:* Because NPStat is not well suited for contigs, and we don't have chromosomes identified in our assembly, we need to concatenate all of our contigs into a single linear "pseudochromosome." This is conceptually the cleanest and simplest approach to follow for NPStat.


- **1.1.** Concatenate the contigs into a single sequence:
   ``` 
   cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/02_IndexRef/ref_full
   
   cat scaffolds_reduced_contigs_kept.fasta | sed -e '1!{/^>.*/d;}' | sed  ':a;N;$!ba;s/\n//2g' > scaffolds_reduced_contigs_kept_concatenated.fasta
   ```

- **1.2.** Check length is correct
  - `02_IndexRef_concat_quast.slurm`
  - Look at the .html output (`02_IndexRef_concat_quast_report.html`) and make sure the concatenated and the original are the same # of bp - they are! (59679359bp)

- **1.3.** Rename the header
  - Just downloaded it and manually redid it and reuploaded it... probably an easier way to do this but nano wasn't working
  - Check that it renamed properly
  ```
  grep -e ">" scaffolds_reduced_contigs_kept_concatenated.fasta
  >scaffolds_reduced_contigs_kept_concatenated
  
  # good to go!
  ```
  
---

### 2. Index the concatenated pseudochromosome reference

- *Concept:* Gota do this for the next step. This is named 'fullconcatref' but it's really just reduced (cleaned) reference that has been concatenated


- **2.1.** Index the reference:
  - `02_IndexRef_fullconcatref.slurm`
  
  ---
  
### 3. Align the reads to the concatenated pseudochromosome reference

- *Concept:* I need to go through SNP steps a, b, and c before considering the .bam files ready to be mpileup'd (by sample). Steps b and c may not be totally necessary since the reads will already be filtered through the -snpfile sites that I point NPStat toward, but I don't think it will hurt anything and will decrease the size of the .mpileup file for NPStat.

- **3.1.** Convert fastq files to bams:
  - `03_SNPa_Fastq2Bam_CERW_reducedrefconcat.slurm`
  - `03_SNPa_Fastq2Bam_PROW_reducedrefconcat.slurm`

- **3.2.** Mark duplicates:
  - `03_SNPb_MrkDup_CERW_reducedrefconcat.slurm`
  - `03_SNPb_MrkDup_PROW_reducedrefconcat.slurm`
  
- **3.3.** Skip alignments with MAPQ smaller than 20, calculate average depth, and index *align_sort_dm_mq20.bams:
  - `03_SNPc_MapQDepth_CERW_reducedrefconcat.slurm`
  - `03_SNPc_MapQDepth_PROW_reducedrefconcat.slurm`
  
  ---

### 4. Generate snpfile

- *Concept:* This -snpfile is to run NPStat on only the SNP variants of interest, not the whole contig/chromosome, etc. However, because the files denoting the variant sites of interest we have currently are based on contig positions (contig | site), they are not necessarily unique and are definitely not the same position as is on the (except for the first contig), so need to calculate a running total while assigning sites to create the single column -snpfile for input into NPStat

- **4.1.** Make .bed/coordinate files that have our sites of interest:
  - Note: I got this information from the PLINK prune.in (intervals) file and did some savvy cp/fr to get these formatted properly: `CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.intervals.forbams_final.bed` and `PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.intervals.forbams_final.bed`

- **4.2.** Prepare snpfile with python scripts
  - Mary DuBose created this python script - **THANK YOU!!!**
    - `generate_NPStat_snpfile_CERW.py`
    - `generate_NPStat_snpfile_PROW.py`
  - It takes in the .bed file out does all the math necessary to make the snpfile with the positions as if they were linear across the entire concatenated reference genome
  - It does not require the node name values be in sequential order (i.e., the .bed file can skip a node number... e.g., node1, node2, node4, node5 will work)
  
  
- **4.3.** Run on AHPCC
  - Adjust the file names/locations within the scripts and then...
  
  ```
  module load python/3.8-anaconda
  python --version # 3.9.6
  cd /scrfs/storage/amatthews/python_scripts
  python generate_NPStat_snpfile_CERW.py
  python generate_NPStat_snpfile_PROW.py
  ```



---

### 5. Generate .pileups by sample

- *Concept:* Need individual .pileup files (individual by sample) to feed into NPStat. I also prepped the samples names files with the ploidy at this step (to reuse in the NPStat loop)

- **5.1.** Prep the sample names/ploidy files:
  - `CERW_ALL_samples_hapsize_NPStat.txt`
  - `PROW_FILTERED_samples_hapsize_NPStat.txt`
  
- **5.2.** Generate the sample pileup files:
  - `04_NPStat_pileupbysample_CERW_reducedrefconcat.slurm`
  - `04_NPStat_pileupbysample_PROW_reducedrefconcat.slurm`
  
  
---

### 6. Run NPStat

- *Concept:* We finally get to do to the good stuff! Luca said to use the genome size minus 1 for the window length flag (-l) because of how we have created this pseudochromosome. Anything less (like 10k or 100k) produces all 0s/NAs - I can confirm this happens!

- *Note on non-loop:* I tried this in a while loop using a slurm and the results were weirdly all 0s/NAs. I'm not sure what could be happening, but because I know it works properly when I do one sample at a time, I'm gonna do one sample at a time!

- **6.1.** Set up each line of code per sample:
  - Note: these two files are what I'm calling the 'true haploidy' value (e.g., 84 mites = 164 ploidy). See `6.4.` for alternate ploidy if needed
  - - `CERW_ALL_samples_hapsize_NPStat.txt`
  - - `PROW_FILTERED_samples_hapsize_NPStat.txt`
  - Use these as guides to make the one-by-one/sample-by-sample code. Double and triple check they are correct!
  
- **6.2.** Run NPStat - true ploidy
  - `04_NPStat_truePloidy_CERW_reducedrefconcat_onebyone.txt` - using the 'true haploidy' value
  - `04_NPStat_truePloidy_PROW_reducedrefconcat_onebyone.txt` - using the 'true haploidy' value
  
- **6.3.** Make results directory and move results into new directory - true ploidy
  ```
  mkdir -p /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_concat_CERW/NPStat_out/true_ploidy-works
  mkdir -p /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_concat_PROW/NPStat_out/true_ploidy-works
  
  # Move true ploidy results over to new dir
  cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_concat_CERW/bam/
  mv *.stats ../NPStat_out/true_ploidy-works
  
  cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_concat_PROW/bam/
  mv *.stats ../NPStat_out/true_ploidy-works
  ```
  
- **6.4.** Get summary file of results - true ploidy
  ```
  # move to proper dir
  cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_concat_CERW/NPStat_out/true_ploidy-works
  
  # Remove first line from all stats files
  sed -i '1d' *.stats
  
  # Cat all results together for easy summary 
  cat *.stats > CERW_NPStat_truePloidy_concat.txt
  
  # move to proper dir
  cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_concat_PROW/NPStat_out/true_ploidy-works
  
  # Remove first line from all stats files
  sed -i '1d' *.stats
  
  # Cat all results together for easy summary 
  cat *.stats > PROW_NPStat_truePloidy_concat.txt
  ```
  
  - Outputs: 
  - - `CERW_NPStat_truePloidy_concat.txt`
  - - `PROW_NPStat_truePloidy_concat.txt`
  - Need to manually remove any window #2s and add header line to these output files and can assign samples (they go in numerical order)



  
