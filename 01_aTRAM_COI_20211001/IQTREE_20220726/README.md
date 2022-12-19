Files used for estimating phylogeny of 57 SNP samples

On my Linux machine:

```
pwd: /home/lse305/Desktop/Alix/snp_20220705/phylogeny

## mafft MAFFT v7.503

mafft --auto SNP-phylo.fasta > SNP-phylo.mafftaligned.fasta

## Strategy: FFT-NS-i (Standard)

## iqtree version 1.6.12

# AICc merit was best
iqtree -s SNP-phylo.mafftaligned.fasta -m MFP -merit AICc -bb 1000 -pre COI_AICc -nt AUTO
```
