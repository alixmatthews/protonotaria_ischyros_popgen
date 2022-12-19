### poolfstat
### Alix Matthews
### 2022 Sept 7

## Running only on ld01, dp15 datasets

### R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"

### OS: Ubuntu 20.04.3 LTS


#### Set working directory and load libraries ####
setwd("~/Desktop/Alix/snp_20220730")

library(poolfstat) # 2.1.1
library(phangorn) # 2.9.0


# following the poolfstat vignette: https://cran.r-project.org/web/packages/poolfstat/vignettes/vignette.pdf
# and https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13557 (supplemental file 3)


#### ~~~ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01  ~~~ ####
# load .vcf and make pooldata object
cerw_dp15_ld01_pfs<-vcf2pooldata("CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf.gz",poolsizes=c(94, 160, 140, 188, 52, 18, 140, 180, 168, 118, 168, 200, 200, 130, 200, 200, 128, 158, 198, 200, 36, 200, 200, 200, 84, 200, 200, 188, 200, 108), poolnames = c("CERW067A_PA","CERW068_PA","CERW069_PA","CERW070_PA","CERW358A_PA","CERW371B_PA","CERW575R_OZ","CERW577R_OZ","CERW578R_OZ","CERW687R_OZ","CERW691R_OZ","CERW702R_OZ","CERW717R_TN","CERW718R_TN","CERW719R_TN","CERW720R_TN","CERW721R_TN","CERW723R_TN","CERW733R_IN","CERW734R_IN","CERW735R_IN","CERW736R_IN","CERW737R_IN","CERW739R_IN","CERW744R_WI","CERW746R_WI","CERW747R_WI","CERW748R_WI","CERW750R_WI","CERW751R_WI"))

# pool size is in haploid seq value... so specifying a pool size of 200 haploid sequences means that pools consisted of 100 diploid individuals

cerw_dp15_ld01_pfs # overview,  * Number of SNPs   =  8735

# Standard-error of  FST estimates can be estimated using a block-jackknife sampling approach (see Appendix A.1) by specifying the number of consecutive SNPs defining a block with the argument nsnp.per.bjack.block (by default nsnp.per.bjack.block=0, i.e., no block-jackknife is carried out) as illustrated below for:

# setting the block to 5 in order to include the most number of SNPs
cerw_dp15_ld01_pfs.fst.jack<-computeFST(cerw_dp15_ld01_pfs,nsnp.per.bjack.block = 10, verbose=TRUE)
# 362 Jackknife blocks identified with 3620 SNPs (out of 8735 ).

cerw_dp15_ld01_pfs.fst.jack<-computeFST(cerw_dp15_ld01_pfs,nsnp.per.bjack.block = 20, verbose=TRUE)
# 44 Jackknife blocks identified with 880 SNPs (out of 8735 ).

cerw_dp15_ld01_pfs.fst.jack<-computeFST(cerw_dp15_ld01_pfs,nsnp.per.bjack.block = 5, verbose=TRUE)
# 1216 Jackknife blocks identified with 6080 SNPs (out of 8735 ).

cerw_dp15_ld01_pfs.fst.jack$FST #genome-wide Fst over all populations: 0.03522641
cerw_dp15_ld01_pfs.fst.jack$mean.fst # #block-jacknife estimate of s.e: 0.03275713 (this value will change based on nsnp.per.bjack.block)
cerw_dp15_ld01_pfs.fst.jack$se.fst #s.e. of the genome-wide Fst estimate: 0.0009387949
cerw_dp15_ld01_pfs.fst.jack$mean.fst+c(-1.96,1.96)*cerw_dp15_ld01_pfs.fst.jack$se.fst #95% c.i. of the estimated genome-wide Fst: 0.03091709 0.03459716. These values, however, do not overlap the genome-wide Fst over all populations as identified above (0.03522641), so our data are the 5% chance that the mean falls outside of the CIs


# 3.2 Estimating and visualizing pairwise-population FST
# 3.2.1 The compute.pairwiseFST and the heatmap functions

cerw_dp15_ld01_pfs.pairwisefst<-compute.pairwiseFST(cerw_dp15_ld01_pfs,verbose=FALSE)
heatmap(cerw_dp15_ld01_pfs.pairwisefst)

# Saved as A4, landscape (cuts off a little bit of the bottom labels, but best I can get it)

# reformatting for import into SplitsTree
cerw_dp15_ld01_pfs.pairwisefst_mat<- as.matrix(cerw_dp15_ld01_pfs.pairwisefst@PairwiseFSTmatrix)
class(cerw_dp15_ld01_pfs.pairwisefst_mat)
cerw_dp15_ld01_pfs.pairwisefst_dist<-as.dist(cerw_dp15_ld01_pfs.pairwisefst_mat)
class(cerw_dp15_ld01_pfs.pairwisefst_dist)
write.nexus.dist(cerw_dp15_ld01_pfs.pairwisefst_dist, file = "./nexus_files/cerw_dp15_ld01_pfs.pairwisefst_dist.nexus")

# 3.2.2 Block-Jackknife estimation of FˆST standard-error and visualisation of confidence intervals

cerw_dp15_ld01_pfs.pairwise.fst.jack<-compute.pairwiseFST(cerw_dp15_ld01_pfs, nsnp.per.bjack.block = 50, verbose=FALSE)
head(cerw_dp15_ld01_pfs.pairwise.fst.jack@values)
plot(cerw_dp15_ld01_pfs.pairwise.fst.jack)
write.csv(cerw_dp15_ld01_pfs.pairwise.fst.jack@PairwiseFSTmatrix, "./pairwise_fst_matrices/cerw_dp15_ld01_pfs.pairwise.fst.jack_pcoa.csv" )


# Pool subsets
cerw_dp15_ld01_pfs.PA.subset<-pooldata.subset(cerw_dp15_ld01_pfs,pool.index=c(1:6),verbose=FALSE)
cerw_dp15_ld01_pfs.PA.subset #display a summary of the resulting pooldata object

cerw_dp15_ld01_pfs.OZ.subset<-pooldata.subset(cerw_dp15_ld01_pfs,pool.index=c(7:12),verbose=FALSE)
cerw_dp15_ld01_pfs.OZ.subset #display a summary of the resulting pooldata object

cerw_dp15_ld01_pfs.TN.subset<-pooldata.subset(cerw_dp15_ld01_pfs,pool.index=c(13:18),verbose=FALSE)
cerw_dp15_ld01_pfs.TN.subset #display a summary of the resulting pooldata object

cerw_dp15_ld01_pfs.IN.subset<-pooldata.subset(cerw_dp15_ld01_pfs,pool.index=c(19:24),verbose=FALSE)
cerw_dp15_ld01_pfs.IN.subset #display a summary of the resulting pooldata object

cerw_dp15_ld01_pfs.WI.subset<-pooldata.subset(cerw_dp15_ld01_pfs,pool.index=c(25:30),verbose=FALSE)
cerw_dp15_ld01_pfs.WI.subset #display a summary of the resulting pooldata object


# genome wide Fst per population
cerw_dp15_ld01_pfs.PA.fst<-computeFST(cerw_dp15_ld01_pfs.PA.subset)
cerw_dp15_ld01_pfs.PA.fst$FST  # 0.07324691

cerw_dp15_ld01_pfs.PA.fst.jack<-computeFST(cerw_dp15_ld01_pfs.PA.subset,nsnp.per.bjack.block = 5, verbose=TRUE)
cerw_dp15_ld01_pfs.PA.fst.jack$FST #0.07324691
cerw_dp15_ld01_pfs.PA.fst.jack$mean.fst #0.06817799
cerw_dp15_ld01_pfs.PA.fst.jack$se.fst #0.001999233

 
cerw_dp15_ld01_pfs.OZ.fst<-computeFST(cerw_dp15_ld01_pfs.OZ.subset)
cerw_dp15_ld01_pfs.OZ.fst$FST  # 0.02233612

cerw_dp15_ld01_pfs.OZ.fst.jack<-computeFST(cerw_dp15_ld01_pfs.OZ.subset,nsnp.per.bjack.block = 5, verbose=TRUE)
cerw_dp15_ld01_pfs.OZ.fst.jack$FST #0.02233612
cerw_dp15_ld01_pfs.OZ.fst.jack$mean.fst #0.02120492
cerw_dp15_ld01_pfs.OZ.fst.jack$se.fst #0.0008320492


cerw_dp15_ld01_pfs.TN.fst<-computeFST(cerw_dp15_ld01_pfs.TN.subset)
cerw_dp15_ld01_pfs.TN.fst$FST  # 0.03738247

cerw_dp15_ld01_pfs.TN.fst.jack<-computeFST(cerw_dp15_ld01_pfs.TN.subset,nsnp.per.bjack.block = 5, verbose=TRUE)
cerw_dp15_ld01_pfs.TN.fst.jack$FST #0.03738247
cerw_dp15_ld01_pfs.TN.fst.jack$mean.fst #0.03434
cerw_dp15_ld01_pfs.TN.fst.jack$se.fst #0.001251473



cerw_dp15_ld01_pfs.IN.fst<-computeFST(cerw_dp15_ld01_pfs.IN.subset)
cerw_dp15_ld01_pfs.IN.fst$FST  # 0.01423587

cerw_dp15_ld01_pfs.IN.fst.jack<-computeFST(cerw_dp15_ld01_pfs.IN.subset,nsnp.per.bjack.block = 5, verbose=TRUE)
cerw_dp15_ld01_pfs.IN.fst.jack$FST #0.01423587
cerw_dp15_ld01_pfs.IN.fst.jack$mean.fst #0.01331207
cerw_dp15_ld01_pfs.IN.fst.jack$se.fst #0.0007756753



cerw_dp15_ld01_pfs.WI.fst<-computeFST(cerw_dp15_ld01_pfs.WI.subset)
cerw_dp15_ld01_pfs.WI.fst$FST  # 0.02326028

cerw_dp15_ld01_pfs.WI.fst.jack<-computeFST(cerw_dp15_ld01_pfs.WI.subset,nsnp.per.bjack.block = 5, verbose=TRUE)
cerw_dp15_ld01_pfs.WI.fst.jack$FST #0.02326028
cerw_dp15_ld01_pfs.WI.fst.jack$mean.fst #0.02171802
cerw_dp15_ld01_pfs.WI.fst.jack$se.fst #0.00105776



# pairwise fst within pops
cerw_dp15_ld01_pfs.PA.pairwisefst<-compute.pairwiseFST(cerw_dp15_ld01_pfs.PA.subset, nsnp.per.bjack.block = 5, verbose=FALSE)
heatmap(cerw_dp15_ld01_pfs.PA.pairwisefst)
plot_fstats(cerw_dp15_ld01_pfs.PA.pairwisefst)

cerw_dp15_ld01_pfs.OZ.pairwisefst<-compute.pairwiseFST(cerw_dp15_ld01_pfs.OZ.subset,nsnp.per.bjack.block = 5,verbose=FALSE)
heatmap(cerw_dp15_ld01_pfs.OZ.pairwisefst)
plot_fstats(cerw_dp15_ld01_pfs.OZ.pairwisefst)

cerw_dp15_ld01_pfs.TN.pairwisefst<-compute.pairwiseFST(cerw_dp15_ld01_pfs.TN.subset,nsnp.per.bjack.block = 5,verbose=FALSE)
heatmap(cerw_dp15_ld01_pfs.TN.pairwisefst)
plot_fstats(cerw_dp15_ld01_pfs.TN.pairwisefst)

cerw_dp15_ld01_pfs.IN.pairwisefst<-compute.pairwiseFST(cerw_dp15_ld01_pfs.IN.subset,nsnp.per.bjack.block = 5,verbose=FALSE)
heatmap(cerw_dp15_ld01_pfs.IN.pairwisefst)
plot_fstats(cerw_dp15_ld01_pfs.IN.pairwisefst)

cerw_dp15_ld01_pfs.WI.pairwisefst<-compute.pairwiseFST(cerw_dp15_ld01_pfs.WI.subset,nsnp.per.bjack.block = 5,verbose=FALSE)
heatmap(cerw_dp15_ld01_pfs.WI.pairwisefst)
plot_fstats(cerw_dp15_ld01_pfs.WI.pairwisefst)


































#### ~~~ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01  ~~~ ####
# load .vcf and make pooldata object
prow_dp15_ld01_pfs<-vcf2pooldata("PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf.gz",poolsizes=c(200, 180, 200, 200, 200, 120, 200, 118, 98, 98, 200, 132, 184, 114, 146, 200, 200, 140, 200, 200, 200, 200, 200, 200, 166, 94, 200), poolnames = c("PROW671R_AR","PROW672R_AR","PROW684R_AR","PROW685R_AR","PROW690R_AR","PROW709R_AR","PROW752R_WI","PROW753R_WI","PROW754R_WI","PROW756R_WI","PROW758R_WI","PROW759R_WI","PROW805R_OH","PROW807R_OH","PROW808R_OH","PROW809R_OH","PROW810R_OH","PROW813R_OH","PROW825R_SC","PROW828R_SC","PROW829R_SC","PROW831R_SC","PROW833R_SC","PROW974R_VA","PROW976R_VA","PROW977R_VA","PROW978R_VA"))

# pool size is in haploid seq value... so specifying a pool size of 200 haploid sequences means that pools consisted of 100 diploid individuals

prow_dp15_ld01_pfs # overview,  * Number of SNPs   =  6713

# Standard-error of  FST estimates can be estimated using a block-jackknife sampling approach (see Appendix A.1) by specifying the number of consecutive SNPs defining a block with the argument nsnp.per.bjack.block (by default nsnp.per.bjack.block=0, i.e., no block-jackknife is carried out) as illustrated below for:


# set blockto 2 to keep the most number of SNPs, but check others just in case. the genome-wide Fst does not change, but the mean.fst and the se and the CI values change
prow_dp15_ld01_pfs.fst.jack<-computeFST(prow_dp15_ld01_pfs,nsnp.per.bjack.block = 20, verbose=TRUE)
#28 Jackknife blocks identified with 560 SNPs (out of 6713 ).

prow_dp15_ld01_pfs.fst.jack<-computeFST(prow_dp15_ld01_pfs,nsnp.per.bjack.block = 10, verbose=TRUE)
#176 Jackknife blocks identified with 1760 SNPs (out of 6713 ).

prow_dp15_ld01_pfs.fst.jack<-computeFST(prow_dp15_ld01_pfs,nsnp.per.bjack.block = 5, verbose=TRUE)
# 781 Jackknife blocks identified with 3905 SNPs (out of 6713 ).

prow_dp15_ld01_pfs.fst.jack<-computeFST(prow_dp15_ld01_pfs,nsnp.per.bjack.block = 2, verbose=TRUE)
# 2871 Jackknife blocks identified with 5742 SNPs (out of 6713 ).

prow_dp15_ld01_pfs.fst.jack$FST #genome-wide Fst over all populations: 0.06259998
prow_dp15_ld01_pfs.fst.jack$mean.fst # mean fst: 0.06033158
prow_dp15_ld01_pfs.fst.jack$se.fst #s.e. of the genome-wide Fst estimate: 0.001254795
prow_dp15_ld01_pfs.fst.jack$mean.fst+c(-1.96,1.96)*prow_dp15_ld01_pfs.fst.jack$se.fst #95% c.i. of the estimated genome-wide Fst: 0.05787218 0.06279098 # this one does cover the true genome-wide Fst mean, but not that it is not always the case (see CERW for example)


# 3.2 Estimating and visualizing pairwise-population FST
# 3.2.1 The compute.pairwiseFST and the heatmap functions

prow_dp15_ld01_pfs.pairwisefst<-compute.pairwiseFST(prow_dp15_ld01_pfs,verbose=FALSE)
heatmap(prow_dp15_ld01_pfs.pairwisefst)

# Saved as A4, landscape (cuts off a little bit of the bottom labels, but best I can get it)

# reformatting for import into SplitsTree
prow_dp15_ld01_pfs.pairwisefst_mat<- as.matrix(prow_dp15_ld01_pfs.pairwisefst@PairwiseFSTmatrix)
class(prow_dp15_ld01_pfs.pairwisefst_mat)
prow_dp15_ld01_pfs.pairwisefst_dist<-as.dist(prow_dp15_ld01_pfs.pairwisefst_mat)
class(prow_dp15_ld01_pfs.pairwisefst_dist)
write.nexus.dist(prow_dp15_ld01_pfs.pairwisefst_dist, file = "./nexus_files/prow_dp15_ld01_pfs.pairwisefst_dist.nexus")

# 3.2.2 Block-Jackknife estimation of FˆST standard-error and visualisation of confidence intervals

prow_dp15_ld01_pfs.pairwise.fst.jack<-compute.pairwiseFST(prow_dp15_ld01_pfs, nsnp.per.bjack.block = 40, verbose=FALSE)

# error (below) if I use 50 for nsnp.per.bjack.block 
## Error in generate.jackknife.blocks(x, nsnp.per.bjack.block, verbose = verbose) : 
## Exit function: No contig available after applying filtering steps (e.g., try lowering nsnp.per.bjack.block)

head(prow_dp15_ld01_pfs.pairwise.fst.jack@values)
plot(prow_dp15_ld01_pfs.pairwise.fst.jack)
write.csv(prow_dp15_ld01_pfs.pairwise.fst.jack@PairwiseFSTmatrix, "./pairwise_fst_matrices/prow_dp15_ld01_pfs.pairwise.fst.jack_pcoa.csv" )


# Pool subsets
prow_dp15_ld01_pfs.AR.subset<-pooldata.subset(prow_dp15_ld01_pfs,pool.index=c(1:6),verbose=FALSE)
prow_dp15_ld01_pfs.AR.subset #display a summary of the resulting pooldata object

prow_dp15_ld01_pfs.WI.subset<-pooldata.subset(prow_dp15_ld01_pfs,pool.index=c(7:12),verbose=FALSE)
prow_dp15_ld01_pfs.WI.subset #display a summary of the resulting pooldata object

prow_dp15_ld01_pfs.OH.subset<-pooldata.subset(prow_dp15_ld01_pfs,pool.index=c(13:18),verbose=FALSE)
prow_dp15_ld01_pfs.OH.subset #display a summary of the resulting pooldata object

prow_dp15_ld01_pfs.SC.subset<-pooldata.subset(prow_dp15_ld01_pfs,pool.index=c(19:23),verbose=FALSE)
prow_dp15_ld01_pfs.SC.subset #display a summary of the resulting pooldata object

prow_dp15_ld01_pfs.VA.subset<-pooldata.subset(prow_dp15_ld01_pfs,pool.index=c(24:27),verbose=FALSE)
prow_dp15_ld01_pfs.VA.subset #display a summary of the resulting pooldata object


# genome wide Fst per population
prow_dp15_ld01_pfs.AR.fst<-computeFST(prow_dp15_ld01_pfs.AR.subset)
prow_dp15_ld01_pfs.AR.fst$FST  # 0.04095089

prow_dp15_ld01_pfs.AR.fst.jack<-computeFST(prow_dp15_ld01_pfs.AR.subset,nsnp.per.bjack.block = 2, verbose=TRUE)
prow_dp15_ld01_pfs.AR.fst.jack$FST #0.04095089
prow_dp15_ld01_pfs.AR.fst.jack$mean.fst #0.03935842
prow_dp15_ld01_pfs.AR.fst.jack$se.fst #0.001184944



prow_dp15_ld01_pfs.WI.fst<-computeFST(prow_dp15_ld01_pfs.WI.subset)
prow_dp15_ld01_pfs.WI.fst$FST  # 0.04324778

prow_dp15_ld01_pfs.WI.fst.jack<-computeFST(prow_dp15_ld01_pfs.WI.subset,nsnp.per.bjack.block = 2, verbose=TRUE)
prow_dp15_ld01_pfs.WI.fst.jack$FST #0.04324778
prow_dp15_ld01_pfs.WI.fst.jack$mean.fst #0.04155359
prow_dp15_ld01_pfs.WI.fst.jack$se.fst #0.001572207



prow_dp15_ld01_pfs.OH.fst<-computeFST(prow_dp15_ld01_pfs.OH.subset)
prow_dp15_ld01_pfs.OH.fst$FST  # 0.1026248

prow_dp15_ld01_pfs.OH.fst.jack<-computeFST(prow_dp15_ld01_pfs.OH.subset,nsnp.per.bjack.block = 2, verbose=TRUE)
prow_dp15_ld01_pfs.OH.fst.jack$FST #0.1026248
prow_dp15_ld01_pfs.OH.fst.jack$mean.fst #0.09896527
prow_dp15_ld01_pfs.OH.fst.jack$se.fst #0.002214007



prow_dp15_ld01_pfs.SC.fst<-computeFST(prow_dp15_ld01_pfs.SC.subset)
prow_dp15_ld01_pfs.SC.fst$FST  # 0.03655796

prow_dp15_ld01_pfs.SC.fst.jack<-computeFST(prow_dp15_ld01_pfs.SC.subset,nsnp.per.bjack.block = 2, verbose=TRUE)
prow_dp15_ld01_pfs.SC.fst.jack$FST #0.03655796
prow_dp15_ld01_pfs.SC.fst.jack$mean.fst #0.03520399
prow_dp15_ld01_pfs.SC.fst.jack$se.fst #0.001335409



prow_dp15_ld01_pfs.VA.fst<-computeFST(prow_dp15_ld01_pfs.VA.subset)
prow_dp15_ld01_pfs.VA.fst$FST  # 0.09029264

prow_dp15_ld01_pfs.VA.fst.jack<-computeFST(prow_dp15_ld01_pfs.VA.subset,nsnp.per.bjack.block = 2, verbose=TRUE)
prow_dp15_ld01_pfs.VA.fst.jack$FST #0.09029264
prow_dp15_ld01_pfs.VA.fst.jack$mean.fst #0.08671451
prow_dp15_ld01_pfs.VA.fst.jack$se.fst #0.002130058




# pairwise fst within pops
prow_dp15_ld01_pfs.AR.pairwisefst<-compute.pairwiseFST(prow_dp15_ld01_pfs.AR.subset,nsnp.per.bjack.block = 2,verbose=FALSE)
heatmap(prow_dp15_ld01_pfs.AR.pairwisefst)
plot_fstats(prow_dp15_ld01_pfs.AR.pairwisefst)

prow_dp15_ld01_pfs.WI.pairwisefst<-compute.pairwiseFST(prow_dp15_ld01_pfs.WI.subset,nsnp.per.bjack.block = 2,verbose=FALSE)
heatmap(prow_dp15_ld01_pfs.WI.pairwisefst)
plot_fstats(prow_dp15_ld01_pfs.WI.pairwisefst)

prow_dp15_ld01_pfs.OH.pairwisefst<-compute.pairwiseFST(prow_dp15_ld01_pfs.OH.subset,nsnp.per.bjack.block = 2,verbose=FALSE)
heatmap(prow_dp15_ld01_pfs.OH.pairwisefst)
plot_fstats(prow_dp15_ld01_pfs.OH.pairwisefst)

prow_dp15_ld01_pfs.SC.pairwisefst<-compute.pairwiseFST(prow_dp15_ld01_pfs.SC.subset,nsnp.per.bjack.block = 2,verbose=FALSE)
heatmap(prow_dp15_ld01_pfs.SC.pairwisefst)
plot_fstats(prow_dp15_ld01_pfs.SC.pairwisefst)

prow_dp15_ld01_pfs.VA.pairwisefst<-compute.pairwiseFST(prow_dp15_ld01_pfs.VA.subset,nsnp.per.bjack.block = 2,verbose=FALSE)
heatmap(prow_dp15_ld01_pfs.VA.pairwisefst)
plot_fstats(prow_dp15_ld01_pfs.VA.pairwisefst)

