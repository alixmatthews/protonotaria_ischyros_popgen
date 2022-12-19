These scripts use aTRAM to retrieve COI from the 30 CERW samples and 30 PROW samples (from SNP/pop gen project)

COI reference is amino acid translated COI from Proctophyllodes miliariae

```01_aTRAM_COI_20211001.slurm``` is the first attempt to retrieve COI from the 60 SNP samples. Something went wrong with the exonerate pipeline for every other sample (strangely). So the ```01_aTRAM_COI_20211001_v2.slurm``` is the same slurm with different input file (only includes the 30 samples that did not get properly exonerated), and only the "stitcher"/exonerate pipeline is active (all else is commented out to hopefully save time...)

---
#### Concatenating all the COI .fasta files

- Downloaded all the stitched (exonerate'd) .fasta files for each sample (COI only in this case) on my own computer ... manually ... could probably make a loop to do this but the nested file structure is a little complicated.
- Cat them together: ```cat *.fasta > CERW_and_PROW_snp.fasta```
- Check there are 60 files: ```grep -c '>' CERW_and_PROW_snp.fasta```
- copy and paste them to the big file ```ParulidMites_2017_2021.fasta``` (which also includes my original COI, *Tyrannidectes*, 60 SNP samples, and 33 Kays samples)

---
#### NOTE: This was cp over from 20210816_phylo repo
---

#### IQTREE_20220726

This section was not cp over from the older repo, but I did add it now. Contains original fasta and mafft-aligned fasta as well as iq-tree outputs and figures



---

