## Downstream analyses on the filtered datasets (from the last 03_SNPk step)

### Multiple analyses...
1. `04_plink_LD_PCA` - LD pruning steps
2. `04_ADMIXTURE_20220709`- ADMIXTURE analysis
3. `04_NPStat_20220824` - NPStat
4. Stuff in R...
- Need to run `After LD-pruning...` code below (for import into R)
  - 4i. `04_R_PCA_DAPC_20220730` - DAPC analyses
  - 4ii. `04_R_IBD_20220907` - isolation by distance
  - 4iii. `04_R_poolfstat_20220907` - poolfstat


---

## After LD-pruning...
### need to select LD-pruned sites and gzip files for import into R
### These are 15X (ld01)

#### *A. ischyros* first
```
module load java/sunjdk_1.8.0 # for gatk and picard
module load gatk/4.2.6.1
module load fastqc/0.11.5 # workaround for picard-tools
module load picard-tools/2.17.10
module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate BCFTools

## Select the variants from prune.in (variants to keep; prune.out are variants that are to be removed because potentially linked). So the next step in this section is to select only the prune.in variants to make the .vcf and the gzip this for import into R


## DP = 15 (15X; ld01)
## Select the variants from prune.in (variants to keep; prune.out are variants that are to be removed because potentially linked). So the next step in this section is to select only the prune.in variants to make the .vcf and the gzip this for import into R

## Duplicating/renaming the prune.in file (variants to keep, file looks like 'CHR:POS') as .intervals for GATK to recognize

cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_CERW/vcf/CERW_plink_LD_PCA_20220706

cp CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.in CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.intervals

## need to index the maf05_dp15_maxmiss100 vcf (vcf above the ld pruning)
gatk IndexFeatureFile --input /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_CERW/vcf/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf

## then can do SelectVariants
gatk SelectVariants --variant /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_CERW/vcf/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --intervals /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_CERW/vcf/CERW_plink_LD_PCA_20220706/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.intervals --output /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_CERW/vcf/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf

## then gzip these .vcf files
cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_CERW/vcf

bcftools view --output-type z --output CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf.gz CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf

```



#### *A. protonotaria* second

```

##  DP = 15 (15X; ld01)

## Duplicating/renaming the prune.in file (variants to keep, file looks like 'CHR:POS') as .intervals for GATK to recognize

cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_PROW/vcf/PROW_plink_LD_PCA_20220706

cp PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.in PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.intervals


## need to index the maf05_dp15_maxmiss100 vcf (vcf above the ld pruning)
gatk IndexFeatureFile --input /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_PROW/vcf/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf

## then can do SelectVariants
gatk SelectVariants --variant /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_PROW/vcf/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100.recode.vcf --intervals /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_PROW/vcf/PROW_plink_LD_PCA_20220706/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.intervals --output /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_PROW/vcf/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf

## then gzip these .vcf files
cd /scrfs/storage/amatthews/20210816_projects/20210816_snp/03_SNP_20220513/RESULTS_ref_reduced_PROW/vcf

bcftools view --output-type z --output PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf.gz PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf

```

