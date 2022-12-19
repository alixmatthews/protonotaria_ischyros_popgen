Pre-processing was conducted over two rounds of sequencing. The same methods were used for both batches.

**Batch 1.** 20200412

Files needed:  `adapters.fa`, `filenames.txt`
- 1. `00_PP_20210427` - run first
- 2. `00_PP_20210507_bb` - run second


**Batch 2.** 20210902

Files needed: `adapters.fa`, `20210816_snp_filenames.txt`
- 1. `00_PP_20210902_a_snp.slurm` - run first
- 2. `00_PP_20210902_b_snp.slurm` - run second
- 3. `00_PP_SNPDATA_multiqc_before_q_trim.slurm` - multiqc report before quality trimming
