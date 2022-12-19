# SNP calling pipeline

**Note:** for all of these slurms, we are using the "reduced full reference with microbes removed" 
(how we got to that point with that reference assembly is outlined in 20210816_projects/20210816_snp/20220513_snp/02_IndexRef/)

---

## 1. Convert .fastq to indexed, align_sort.bam files 

- slurms: `03_SNPa_Fastq2Bam_CERW_reducedref.slurm` and `03_SNPa_Fastq2Bam_PROW_reducedref.slurm`
- Get flagstat summaries (output: `flagstat_summaries_reducedref.xlsx`)

```
grep "" *.txt > $SPP_flagstat_out_all_files.txt
# import to Excel
# copy and paste first couple of lines and do flash fill to separate #s following : from the sample name
# use filters to filter our what category want to look at and summarize the output
```



## 2. Mark duplicates and index the align_sort_dm.bams
- slurms: `03_SNPb_MrkDup_CERW_reducedref.slurm` and `03_SNPb_MrkDup_PROW_reducedref.slurm`



## 3. Skip alignments with MAPQ smaller than 20, calculate average depth, and index align_sort_dm_mq20.bams 
- slurms: `03_SNPc_MapQDepth_CERW_reducedref.slurm` and `03_SNPc_MapQDepth_PROW_reducedref.slurm`
- summary is here (produced using code below): `03_SNPc_MapQDepth_PROW_CERW_output.xlsx`
- get a list of bam files:


```
find "$PWD" -type f -name '*mq20.bam' > $SPP_align_sort_dm_mq20_bam.list

# grep the avg and stdev: 
grep "Average" *mq20_depth.txt
grep "Stdev" *mq20_depth.txt
```
    


## 4. Calculate depth to get file to import into R and plot 
- slurms: `03_SNPd_DepthOnly_CERW_reducedref.slurm` and `03_SNPd_DepthOnly_PROW_reducedref.slurm`
- NOTE: was not able to import the depth file into R, but the depth (output) files are available to try again



## 5. Create mpileup files, max read depth per file is set at 50
- slurms: `03_SNPe_mpileup_CERW_reducedref.slurm` and `03_SNPe_mpileup_PROW_reducedref.slurm`



## 6. Call variants on the mpileup file, joint-genotyping and rename headers
- slurms: `03_SNPf_callvars_CERW_reducedref.slurm` and `03_SNPf_callvars_PROW_reducedref.slurm`
    + need `CERW_newnames_locations.txt` and `PROW_newnames_locations.txt` with these slurms
- Check the order of the names of each .vcf and then will want to rename and make sure they go in the same order as the original
 ```
 bcftools query -l $SPP_ALL.vcf
 bcftools query -l $SPP_ALL_renamed.vcf
 ```



## 7. Filter by phred QUAL (threshold is 30)
- slurms: `03_SNPg_filterqual_CERW_reducedref.slurm` and `03_SNPg_filterqual_PROW_reducedref.slurm`



## 8. Filter out SNPs that differ from the reference but are equal to each other across samples
- slurms: `03_SNPh_filterminac_CERW_reducedref.slurm` and `03_SNPh_filterminac_PROW_reducedref.slurm`



## 9. Filter on minor allele freq (MAF), several thresholds
- slurms: `03_SNPi_filterminaf_CERW_reducedref.slurm` and `03_SNPi_filterminaf_PROW_reducedref.slurm`



## 10. Filter by sampling seq depth (FORMAT/DP), several thresholds
- slurms: `03_SNPj_filterdp_CERW_reducedref.slurm` and `03_SNPj_filterdp_PROW_reducedref.slurm`


## 11. Remove samples with very low mapping % and poor GC content

```
module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate BCFTools

cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_PROW/vcf

bcftools view --samples ^PROW830R_SC,PROW975R_VA,PROW979R_VA --output PROW_FILTERED_renamed_q30_minac1_maf05_dp05.vcf.gz PROW_ALL_renamed_q30_minac1_maf05_dp05.vcf.gz
bcftools query -l PROW_FILTERED_renamed_q30_minac1_maf05_dp05.vcf.gz

bcftools view --samples ^PROW830R_SC,PROW975R_VA,PROW979R_VA --output PROW_FILTERED_renamed_q30_minac1_maf05_dp10.vcf.gz PROW_ALL_renamed_q30_minac1_maf05_dp10.vcf.gz
bcftools query -l PROW_FILTERED_renamed_q30_minac1_maf05_dp10.vcf.gz

bcftools view --samples ^PROW830R_SC,PROW975R_VA,PROW979R_VA --output PROW_FILTERED_renamed_q30_minac1_maf05_dp15.vcf.gz PROW_ALL_renamed_q30_minac1_maf05_dp15.vcf.gz
bcftools query -l PROW_FILTERED_renamed_q30_minac1_maf05_dp15.vcf.gz

bcftools view --samples ^PROW830R_SC,PROW975R_VA,PROW979R_VA --output PROW_FILTERED_renamed_q30_minac1_maf05_dp20.vcf.gz PROW_ALL_renamed_q30_minac1_maf05_dp20.vcf.gz
bcftools query -l PROW_FILTERED_renamed_q30_minac1_maf05_dp20.vcf.gz

```


## 12. Filter out SNP sites with missing data, several thresholds
- slurms: `03_SNPk_filtermissing_CERW_reducedref.slurm` and `03_SNPk_filtermissing_PROW_reducedref.slurm`


### File with overview of each filtering step: `NumberOfSNPs_filtering.xlsx`


## 13. LD-pruning; move to next directory (`04_Analyses`)
- these slurms are in `04_Analyses/04_plink_LD_PCA_20220617`
