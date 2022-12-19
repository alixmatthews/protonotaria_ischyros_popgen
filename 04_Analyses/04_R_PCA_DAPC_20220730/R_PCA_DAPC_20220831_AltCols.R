### DAPC, PCA, etc. 
### Alix Matthews
### 2022 August 31
##### Same code as R_PCA_DAPC_20220730.R, just different colors and loading using .rda files (not re-running script!)

### R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"

### OS: Ubuntu 20.04.3 LTS

#### Set working directory and load libraries ####
setwd("~/Desktop/Alix/snp_20220730/")


library(vcfR) # 1.12.0 
library(ggplot2) # 3.3.6
library(ape) # 5.6-2
library(RColorBrewer) # 1.1-3
library(ggthemes) # 4.2.4
library(shiny) # 1.7.1
library(adegenet) # 2.1.7
library(poppr) # 2.9.3
library(dplyr) # 1.0.9
library(tidyr) # 1.2.0
library(reshape2) # 1.4.4

#### ~~~ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01  ~~~ ####
#### + housekeeping ####

cerw_dp15_ld01 <- read.vcfR("CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf.gz")

# variant count: 9051

head(cerw_dp15_ld01)
cerw_dp15_ld01_gl <- vcfR2genlight(cerw_dp15_ld01)

# Warning message:
#   In vcfR2genlight(cerw_dp15_ld01) :
#   Found 316 loci with more than two alleles.
# Objects of class genlight only support loci with two alleles.
# 316 loci will be omitted from the genlight object.

# Check some general stats about the genlight object:
ploidy(cerw_dp15_ld01_gl)
cerw_dp15_ld01_gl@gen # this shows there are 8,735 SNP sites
cerw_dp15_ld01_gl@ind.names # check the order of the ind.names

# now assign populations based on the order of the ind.names
pop(cerw_dp15_ld01_gl) <- as.factor(c(rep("Pennsylvania",6), rep("Ozarks",6), rep("Tennessee", 6), rep("Indiana",6), rep("Wisconsin",6)))

cerw_dp15_ld01_gl@pop # check that pop assignment is correct


#### ++ PCA - biallelic only ####
# cerw_dp15_ld01_pca_null <- glPca(cerw_dp15_ld01_gl, nf = NULL) # this will display a screeplot of eigenvalues and user asked for a number of PC axes to retain

cerw_dp15_ld01_pca <- glPca(cerw_dp15_ld01_gl, nf = 3) # select nf (number of principal components to be retained)

# Plot the screeplot
barplot(100*cerw_dp15_ld01_pca$eig/sum(cerw_dp15_ld01_pca$eig), col = heat.colors(50), main="PCA Eigenvalues (CERW, 15X, LD01)")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

# Assign the PCA scores (and populations) for each PC
cerw_dp15_ld01_pca_scores <- as.data.frame(cerw_dp15_ld01_pca$scores)
cerw_dp15_ld01_pca_scores$pop <- pop(cerw_dp15_ld01_gl)

cols <- brewer.pal(n = nPop(cerw_dp15_ld01_gl), name = "Dark2")

cerw_dp15_ld01_pca_null_p <- ggplot(cerw_dp15_ld01_pca_scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2) +
  #geom_text(aes(label=cerw_dp15_ld01_gl$ind.names)) +
  stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ggtitle("CERW PCA with 3 PCs - 15X, LD r2 = 0.1")
cerw_dp15_ld01_pca_null_p


# Warning message:
#   In MASS::cov.trob(data[, vars]) : Probable convergence failure


#### ++ DAPC - biallelic only, 30x cross-validation loop ####

### Looping through 100 mini-xvals (30 xvals) to find the best range to use for the 1000 xvals
### I tried setting the seed to reproduce these results to no avail. I tried all sorts of combinations of set.seed (inside and outside the loop), adding set.seed(NULL) before xval (both inside and outside the loop), setting the seed to a specific number (within the loop - changing each time - and outside the loop), and no results are reproducing even if I do the same thing twice. Even if I do the entire code outside of the loop with the same seed, the results change (even if I set.seed(NULL) and the reset the seed to a specific number). Everything seems to be converging at about the same, at least with my test runs, so I think 100 runs should be good (even if not the exact reproducible results), it should even out...


set.seed(NULL)

cerw_dp15_ld01_dapc_cv_loopi_n.pca<-list()
cerw_dp15_ld01_dapc_cv_loopi_n.da<-list()

for (i in 1:100)
{
  print(i)
  cerw_dp15_ld01_dapc_cv_i <- xvalDapc(tab(cerw_dp15_ld01_gl, NA.method = "mean"), grp = pop(cerw_dp15_ld01_gl), n.rep = 30, parallel = "multicore", ncpus = 6L)
  
  cerw_dp15_ld01_dapc_cv_loopi_n.pca[[i]] <- cerw_dp15_ld01_dapc_cv_i$DAPC$n.pca
  cerw_dp15_ld01_dapc_cv_loopi_n.da[[i]] <- cerw_dp15_ld01_dapc_cv_i$DAPC$n.da
}



# check some basic stats

max(unlist(cerw_dp15_ld01_dapc_cv_loopi_n.pca))
min(unlist(cerw_dp15_ld01_dapc_cv_loopi_n.pca))
mean(unlist(cerw_dp15_ld01_dapc_cv_loopi_n.pca))
median(unlist(cerw_dp15_ld01_dapc_cv_loopi_n.pca))
hist(unlist(cerw_dp15_ld01_dapc_cv_loopi_n.pca))

# only run 'save' if needed to replace previous files
# save(cerw_dp15_ld01_dapc_cv_loopi_n.pca, file='./rda_files/cerw_dp15_ld01_dapc_cv_loopi_n.pca.rda')
# save(cerw_dp15_ld01_dapc_cv_loopi_n.da, file='./rda_files/cerw_dp15_ld01_dapc_cv_loopi_n.da.rda')



# reformat for ggplot histogram

cerw_dp15_ld01_dapc_cv_loopi_n.pca_df<-as.data.frame(matrix(unlist(cerw_dp15_ld01_dapc_cv_loopi_n.pca),nrow=length(cerw_dp15_ld01_dapc_cv_loopi_n.pca),byrow=TRUE))
colnames(cerw_dp15_ld01_dapc_cv_loopi_n.pca_df)[colnames(cerw_dp15_ld01_dapc_cv_loopi_n.pca_df) == 'V1'] <- 'Number_of_PCs'

# only run 'save' if needed to replace previous file
# save(cerw_dp15_ld01_dapc_cv_loopi_n.pca_df, file='./rda_files/cerw_dp15_ld01_dapc_cv_loopi_n.pca_df.rda')


ggplot(cerw_dp15_ld01_dapc_cv_loopi_n.pca_df, aes(x = Number_of_PCs)) +
  geom_bar() +
  scale_x_continuous(breaks=seq(2,25, by = 2)) +
  scale_y_continuous(breaks=seq(0,100, by = 5)) +
  ggtitle("Distribution of PCs after 100 cross-validations (30 reps) \nCERW DP15 LD01") +
  xlab("Number of PCs retained") +
  ylab("Count")+
  theme_bw()

# Saved as A4 landscape 

# Looks like 16 is obviously our best bet, but let's try with more x-vals around the bulk of the distribution (15-17)






#### +++ 1000 xvals centered around # of PCs determined above ####

cerw_dp15_ld01_dapc_cv1000 <- xvalDapc(tab(cerw_dp15_ld01_gl, NA.method = "mean"), grp = pop(cerw_dp15_ld01_gl), n.pca = 15:17, n.rep = 1000, parallel = "multicore", ncpus = 6L)

cerw_dp15_ld01_dapc_cv1000[-1] # 16 PCs

# $n.pca: 16 first PCs of PCA used
# $n.da: 4 discriminant functions saved
# $var (proportion of conserved variance): 0.59

# only run 'save' if needed to replace previous file
# save(cerw_dp15_ld01_dapc_cv1000, file='./rda_files/cerw_dp15_ld01_dapc_cv1000.rda')
load(file='./rda_files/cerw_dp15_ld01_dapc_cv1000.rda')


# This automatically puts the number of PCs determined in the xval
# legends
scatter.dapc(cerw_dp15_ld01_dapc_cv1000$DAPC, cex = 2, legend = TRUE,
             clabel = FALSE, cleg=1.5, posi.leg = "topright", scree.pca = TRUE,
             posi.pca = "topleft", posi.da="none", xax = 1, yax = 2, inset.solid = 1, col=c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), pch=19, bg.inset = "gray97")

# saved as A4, landscape orientation

# no legends
scatter.dapc(cerw_dp15_ld01_dapc_cv1000$DAPC, cex = 2, legend = FALSE,
             clabel = FALSE, cleg=1.5, scree.pca = FALSE,
             posi.da="none", xax = 1, yax = 2, inset.solid = 1, col=c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), pch=19, bg.inset = "gray97")

# saved as A4, landscape orientation

# assignment plot, shows the accuracy of the model to recover the each sample's sample location... cell color indicates the probability of assigning a given sample to the corresponding location (i.e., membership probability). Red = probability of 1 to white = probability of 0. Blue X's indicate original assigned location cluster.
assignplot(cerw_dp15_ld01_dapc_cv1000$DAPC, cex.lab=0.5, pch=4)

# saved as A4, landscape orientation


# here you can choose # of PCs to retain and use the n.da determined from above to make compoplot
cerw_dp15_ld01_dapc_16 <- dapc(cerw_dp15_ld01_gl, n.pca = 16, n.da = 4)

compoplot(cerw_dp15_ld01_dapc_16, col = cols, posi = 'top')
cerw_dp15_ld01_dapc_16.results <- as.data.frame(cerw_dp15_ld01_dapc_16$posterior)
cerw_dp15_ld01_dapc_16.results$pop <- pop(cerw_dp15_ld01_gl)
cerw_dp15_ld01_dapc_16.results$indNames <- rownames(cerw_dp15_ld01_dapc_16.results)
cerw_dp15_ld01_dapc_16.results <- pivot_longer(cerw_dp15_ld01_dapc_16.results, -c(pop, indNames))
colnames(cerw_dp15_ld01_dapc_16.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

cerw_dp15_ld01_dapc_16_p <- ggplot(cerw_dp15_ld01_dapc_16.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') +
  theme_few() +
  scale_fill_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), name = "Assigned\nPopulation") +
  facet_grid(~Original_Pop, scales = "free") +
  labs(y = "Posterior Membership Probability") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  theme(axis.title.y = element_text(face = "bold", size = 16)) +
  theme(axis.text.y = element_text(face = "bold", size = 10)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(strip.text.x = element_text(face = "bold", size = 16)) +
  theme(legend.title = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(face = "bold", size = 14)) +
  theme(legend.position = "bottom")
cerw_dp15_ld01_dapc_16_p

# saved as A4, landscape orientation













#### ++ FIND.CLUSTERS - BIC - loop through K's (1-10) ####

## first check to see how many PCs will need to be retained in the loop (how many to get to 95% of cumulative variance?). Just using diffNgroup so that it doesn't ask me to select the number of clusters. Although, in this case, the number is selects is 2 (even though BIC K = 1 is the lowest)
cerw_dp15_ld01_gl_pcs<- find.clusters(cerw_dp15_ld01_gl, choose.n.clust = FALSE, criterion = "diffNgroup")

## looks like 27 PCs. Select 27 PCs and 2 clusters (will not be using this variable for anything other than checking this)
## add 27 in the n.pca loop below

maxK <- 10
cerw_dp15_ld01_mat <- matrix(nrow=10, ncol=maxK)
colnames(cerw_dp15_ld01_mat) <- 1:ncol(cerw_dp15_ld01_mat)
for(i in 1:nrow(cerw_dp15_ld01_mat)) {
  cerw_dp15_ld01_gl_grp<-find.clusters(cerw_dp15_ld01_gl, n.pca = 27, choose.n.clust = FALSE, max.n.clust = maxK)
  cerw_dp15_ld01_mat[i,] <- cerw_dp15_ld01_gl_grp$Kstat
}

cerw_dp15_ld01_mat_df <- melt(cerw_dp15_ld01_mat)
colnames(cerw_dp15_ld01_mat_df)[1:3] <- c("Group", "K", "BIC")
cerw_dp15_ld01_mat_df$K <- as.factor(cerw_dp15_ld01_mat_df$K)
head(cerw_dp15_ld01_mat_df)

ggplot(cerw_dp15_ld01_mat_df, aes(x=K, y=BIC)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("A. ischyros DP15, ld01") +
  xlab("Number of groups (K)")

# saved as A4, landscape orientation


# The best K = 1, but have to choose 2 for DAPC analysis or else it will fail.

# Do the same thing again but choose K = 2 (n.clust) and check the assignment plots. Keep n.pca as the same as the previous.

cerw_dp15_ld01_gl_grp<-find.clusters(cerw_dp15_ld01_gl, n.pca = 27, n.clust = 2)

# See the group assignments of each sampling location (i.e.,'pop' assigned upon loading in data)
table(pop(cerw_dp15_ld01_gl), cerw_dp15_ld01_gl_grp$grp)

# Compare K = 2 assignment with pop

table.value(table(pop(cerw_dp15_ld01_gl), cerw_dp15_ld01_gl_grp$grp), col.lab=paste("inferred", 1:2),
            row.lab=paste("original", 1:5))




#### ++ FIND.CLUSTERS - DAPC - biallelic only, 30x cross-validation loop ####
set.seed(NULL)

cerw_dp15_ld01_finddapc_cv_loopi_n.pca<-list()
cerw_dp15_ld01_finddapc_cv_loopi_n.da<-list()

for (i in 1:100)
{
  print(i)
  cerw_dp15_ld01_finddapc_cv_i <- xvalDapc(tab(cerw_dp15_ld01_gl, NA.method = "mean"), grp = cerw_dp15_ld01_gl_grp$grp, n.rep = 30, parallel = "multicore", ncpus = 6L)
  
  cerw_dp15_ld01_finddapc_cv_loopi_n.pca[[i]] <- cerw_dp15_ld01_finddapc_cv_i$DAPC$n.pca
  cerw_dp15_ld01_finddapc_cv_loopi_n.da[[i]] <- cerw_dp15_ld01_finddapc_cv_i$DAPC$n.da
}


# check some basic stats

max(unlist(cerw_dp15_ld01_finddapc_cv_loopi_n.pca))
min(unlist(cerw_dp15_ld01_finddapc_cv_loopi_n.pca))
mean(unlist(cerw_dp15_ld01_finddapc_cv_loopi_n.pca))
median(unlist(cerw_dp15_ld01_finddapc_cv_loopi_n.pca))
hist(unlist(cerw_dp15_ld01_finddapc_cv_loopi_n.pca))

# only run 'save' if needed to replace previous files
# save(cerw_dp15_ld01_finddapc_cv_loopi_n.pca, file='./rda_files/cerw_dp15_ld01_finddapc_cv_loopi_n.pca.rda')
# save(cerw_dp15_ld01_finddapc_cv_loopi_n.da, file='./rda_files/cerw_dp15_ld01_finddapc_cv_loopi_n.da.rda')



# reformat for ggplot histogram

cerw_dp15_ld01_finddapc_cv_loopi_n.pca_df<-as.data.frame(matrix(unlist(cerw_dp15_ld01_finddapc_cv_loopi_n.pca),nrow=length(cerw_dp15_ld01_finddapc_cv_loopi_n.pca),byrow=TRUE))
colnames(cerw_dp15_ld01_finddapc_cv_loopi_n.pca_df)[colnames(cerw_dp15_ld01_finddapc_cv_loopi_n.pca_df) == 'V1'] <- 'Number_of_PCs'

# only run 'save' if needed to replace previous file
# save(cerw_dp15_ld01_finddapc_cv_loopi_n.pca_df, file='./rda_files/cerw_dp15_ld01_finddapc_cv_loopi_n.pca_df.rda')


ggplot(cerw_dp15_ld01_finddapc_cv_loopi_n.pca_df, aes(x = Number_of_PCs)) +
  geom_bar() +
  scale_x_continuous(breaks=seq(2,25, by = 2)) +
  scale_y_continuous(breaks=seq(0,100, by = 5)) +
  ggtitle("Distribution of PCs after 100 cross-validations (30 reps) \nCERW DP15 LD01 using find.clusters") +
  xlab("Number of PCs retained") +
  ylab("Count")+
  theme_bw()

# Saved as A4 landscape 

# Looks like 1 is obviously our best bet, but let's try with more x-vals around the bulk of the distribution (1-2)






#### +++ 1000 xvals centered around # of PCs determined above ####

cerw_dp15_ld01_finddapc_cv1000 <- xvalDapc(tab(cerw_dp15_ld01_gl, NA.method = "mean"), grp = cerw_dp15_ld01_gl_grp$grp, n.pca = 1:2, n.rep = 1000, parallel = "multicore", ncpus = 6L)

cerw_dp15_ld01_finddapc_cv1000[-1] # xx PCs

# $n.pca: 1 first PCs of PCA used
# $n.da: 1 discriminant functions saved
# $var (proportion of conserved variance): 0.046

# only run 'save' if needed to replace previous file
# save(cerw_dp15_ld01_finddapc_cv1000, file='./rda_files/cerw_dp15_ld01_finddapc_cv1000.rda')
load(file='./rda_files/cerw_dp15_ld01_finddapc_cv1000.rda')


# This automatically puts the number of PCs determined in the xval
# legends
scatter.dapc(cerw_dp15_ld01_finddapc_cv1000$DAPC, cex = 2, legend = TRUE,
             clabel = FALSE, cleg=1.5, posi.leg = "topright", scree.pca = TRUE,
             posi.pca = "topleft", posi.da="none", xax = 1, yax = 2, inset.solid = 1, col=c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), pch=19, bg.inset = "gray97")

# saved as A4, landscape orientation

# no legends
scatter.dapc(cerw_dp15_ld01_finddapc_cv1000$DAPC, cex = 2, legend = FALSE,
             clabel = FALSE, cleg=1.5, scree.pca = FALSE,
             posi.da="none", xax = 1, yax = 2, inset.solid = 1, col=c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), pch=19, bg.inset = "gray97")

# saved as A4, landscape orientation

# assignment plot, shows the accuracy of the model to recover the each sample's sample location... cell color indicates the probability of assigning a given sample to the corresponding location (i.e., membership probability). Red = probability of 1 to white = probability of 0. Blue X's indicate original assigned location cluster.
assignplot(cerw_dp15_ld01_finddapc_cv1000$DAPC, cex.lab=0.5, pch=4)

# saved as A4, landscape orientation


# here you can choose # of PCs to retain and use the n.da determined from above to make compoplot
cerw_dp15_ld01_finddapc_1 <- dapc(cerw_dp15_ld01_gl, n.pca = 1, n.da = 1)

compoplot(cerw_dp15_ld01_finddapc_1, col = cols, posi = 'top')
cerw_dp15_ld01_finddapc_1.results <- as.data.frame(cerw_dp15_ld01_finddapc_1$posterior)
cerw_dp15_ld01_finddapc_1.results$pop <- pop(cerw_dp15_ld01_gl)
cerw_dp15_ld01_finddapc_1.results$indNames <- rownames(cerw_dp15_ld01_finddapc_1.results)
cerw_dp15_ld01_finddapc_1.results <- pivot_longer(cerw_dp15_ld01_finddapc_1.results, -c(pop, indNames))
colnames(cerw_dp15_ld01_finddapc_1.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

cerw_dp15_ld01_finddapc_1_p <- ggplot(cerw_dp15_ld01_finddapc_1.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') +
  theme_few() +
  scale_fill_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), name = "Assigned\nPopulation") +
  facet_grid(~Original_Pop, scales = "free") +
  labs(y = "Posterior Membership Probability") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  theme(axis.title.y = element_text(face = "bold", size = 16)) +
  theme(axis.text.y = element_text(face = "bold", size = 10)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(strip.text.x = element_text(face = "bold", size = 16)) +
  theme(legend.title = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(face = "bold", size = 14)) +
  theme(legend.position = "bottom")
cerw_dp15_ld01_finddapc_1_p

# saved as A4, landscape orientation
















































#### ~~~ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01 ~~~ ####
#### + housekeeping ####

prow_dp15_ld01 <- read.vcfR("PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.prune.vcf.gz")

# variant count: 6947

head(prow_dp15_ld01)
prow_dp15_ld01_gl <- vcfR2genlight(prow_dp15_ld01)

# Warning message:
#   In vcfR2genlight(prow_dp15_ld01) :
#   Found 234 loci with more than two alleles.
# Objects of class genlight only support loci with two alleles.
# 234 loci will be omitted from the genlight object.

# Check some general stats about the genlight object:
ploidy(prow_dp15_ld01_gl)
prow_dp15_ld01_gl@gen # this shows there are 6,713 SNP sites
prow_dp15_ld01_gl@ind.names # check the order of the ind.names

# now assign populations based on the order of the ind.names
pop(prow_dp15_ld01_gl) <- as.factor(c(rep("Arkansas",6), rep("Wisconsin",6), rep("Ohio", 6), rep("South Carolina",5), rep("Virginia",4)))

prow_dp15_ld01_gl@pop # check that pop assignment is correct


#### ++ PCA - biallelic only ####
# prow_dp15_ld01_pca_null <- glPca(prow_dp15_ld01_gl, nf = NULL) # this will display a screeplot of eigenvalues and user asked for a number of PC axes to retain

prow_dp15_ld01_pca <- glPca(prow_dp15_ld01_gl, nf = 3) # select nf (number of principal components to be retained)

# Plot the screeplot
barplot(100*prow_dp15_ld01_pca$eig/sum(prow_dp15_ld01_pca$eig), col = heat.colors(50), main="PCA Eigenvalues (PROW)")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

# Assign the PCA scores (and populations) for each PC
prow_dp15_ld01_pca_scores <- as.data.frame(prow_dp15_ld01_pca$scores)
prow_dp15_ld01_pca_scores$pop <- pop(prow_dp15_ld01_gl)

cols <- brewer.pal(n = nPop(prow_dp15_ld01_gl), name = "Dark2")


prow_dp15_ld01_pca_null_p <- ggplot(prow_dp15_ld01_pca_scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2) +
  #geom_text(aes(label=prow_dp15_ld01_gl$ind.names)) +
  stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ggtitle("PROW PCA with 3 PCs - DP=15, LD r2 = 0.1")
prow_dp15_ld01_pca_null_p




#### ++ DAPC - biallelic only, 30x cross-validation loop ####
### Looping through 100 mini-xvals (30 xvals) to find the best range to use for the 1000 xvals

set.seed(NULL)

prow_dp15_ld01_dapc_cv_loopi_n.pca<-list()
prow_dp15_ld01_dapc_cv_loopi_n.da<-list()

for (i in 1:100)
{
  print(i)
  prow_dp15_ld01_dapc_cv_i <- xvalDapc(tab(prow_dp15_ld01_gl, NA.method = "mean"), grp = pop(prow_dp15_ld01_gl), n.rep = 30, parallel = "multicore", ncpus = 6L)
  
  prow_dp15_ld01_dapc_cv_loopi_n.pca[[i]] <- prow_dp15_ld01_dapc_cv_i$DAPC$n.pca
  prow_dp15_ld01_dapc_cv_loopi_n.da[[i]] <- prow_dp15_ld01_dapc_cv_i$DAPC$n.da
}

# check some basic stats

max(unlist(prow_dp15_ld01_dapc_cv_loopi_n.pca))
min(unlist(prow_dp15_ld01_dapc_cv_loopi_n.pca))
mean(unlist(prow_dp15_ld01_dapc_cv_loopi_n.pca))
median(unlist(prow_dp15_ld01_dapc_cv_loopi_n.pca))
hist(unlist(prow_dp15_ld01_dapc_cv_loopi_n.pca))

# only run 'save' if needed to replace previous files
# save(prow_dp15_ld01_dapc_cv_loopi_n.pca, file='./rda_files/prow_dp15_ld01_dapc_cv_loopi_n.pca.rda')
# save(prow_dp15_ld01_dapc_cv_loopi_n.da, file='./rda_files/prow_dp15_ld01_dapc_cv_loopi_n.da.rda')



# reformat for ggplot histogram

prow_dp15_ld01_dapc_cv_loopi_n.pca_df<-as.data.frame(matrix(unlist(prow_dp15_ld01_dapc_cv_loopi_n.pca),nrow=length(prow_dp15_ld01_dapc_cv_loopi_n.pca),byrow=TRUE))
colnames(prow_dp15_ld01_dapc_cv_loopi_n.pca_df)[colnames(prow_dp15_ld01_dapc_cv_loopi_n.pca_df) == 'V1'] <- 'Number_of_PCs'

# only run 'save' if needed to replace previous file
# save(prow_dp15_ld01_dapc_cv_loopi_n.pca_df, file='./rda_files/prow_dp15_ld01_dapc_cv_loopi_n.pca_df.rda')


ggplot(prow_dp15_ld01_dapc_cv_loopi_n.pca_df, aes(x = Number_of_PCs)) +
  geom_bar() +
  scale_x_continuous(breaks=seq(2,25, by = 2)) +
  scale_y_continuous(breaks=seq(0,100, by = 5)) +
  ggtitle("Distribution of PCs after 100 cross-validations (30 reps) \nPROW DP15 LD01") +
  xlab("Number of PCs retained") +
  ylab("Count")+
  theme_bw()

# Looks like 10 is obviously our best bet, but let's try with more x-vals around the bulk of the distribution (9-11)





#### +++ 1000 xvals centered around # of PCs determined above ####

prow_dp15_ld01_dapc_cv1000 <- xvalDapc(tab(prow_dp15_ld01_gl, NA.method = "mean"), grp = pop(prow_dp15_ld01_gl), n.pca = 9:11, n.rep = 1000, parallel = "multicore", ncpus = 6L)

prow_dp15_ld01_dapc_cv1000[-1] # 9 PCs

# $n.pca: 10 first PCs of PCA used
# $n.da: 4 discriminant functions saved
# $var (proportion of conserved variance): 0.427

# only run 'save' if needed to replace previous file
# save(prow_dp15_ld01_dapc_cv1000, file='./rda_files/prow_dp15_ld01_dapc_cv1000.rda')
load(file='./rda_files/prow_dp15_ld01_dapc_cv1000.rda')


# This automatically puts the number of PCs determined in the xval
# legends
scatter.dapc(prow_dp15_ld01_dapc_cv1000$DAPC, cex = 2, legend = TRUE,
             clabel = FALSE, cleg=1.5, posi.leg = "topright", scree.pca = TRUE,
             posi.pca = "topleft", posi.da="none", xax = 1, yax = 2, inset.solid = 1, col=c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), pch=19, bg.inset = "gray97")

# saved as A4, landscape orientation

# no legends
scatter.dapc(prow_dp15_ld01_dapc_cv1000$DAPC, cex = 2, legend = FALSE,
             clabel = FALSE, cleg=1.5, scree.pca = FALSE,
             posi.da="none", xax = 1, yax = 2, inset.solid = 1, col=c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), pch=19, bg.inset = "gray97")

# saved as A4, landscape orientation

# assignment plot, shows the accuracy of the model to recover the each sample's sample location... cell color indicates the probability of assigning a given sample to the corresponding location (i.e., membership probability). Red = probability of 1 to white = probability of 0. Blue X's indicate original assigned location cluster.
assignplot(prow_dp15_ld01_dapc_cv1000$DAPC, cex.lab=0.5, pch=4)

# saved as A4, landscape orientation


# here you can choose # of PCs to retain and use the n.da determined from above to make compoplot
prow_dp15_ld01_dapc_10 <- dapc(prow_dp15_ld01_gl, n.pca = 10, n.da = 4)

compoplot(prow_dp15_ld01_dapc_10, col = cols, posi = 'top')
prow_dp15_ld01_dapc_10.results <- as.data.frame(prow_dp15_ld01_dapc_10$posterior)
prow_dp15_ld01_dapc_10.results$pop <- pop(prow_dp15_ld01_gl)
prow_dp15_ld01_dapc_10.results$indNames <- rownames(prow_dp15_ld01_dapc_10.results)
prow_dp15_ld01_dapc_10.results <- pivot_longer(prow_dp15_ld01_dapc_10.results, -c(pop, indNames))
colnames(prow_dp15_ld01_dapc_10.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

prow_dp15_ld01_dapc_10_p <- ggplot(prow_dp15_ld01_dapc_10.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') +
  theme_few() +
  scale_fill_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), name = "Assigned\nPopulation") +
  facet_grid(~Original_Pop, scales = "free") +
  labs(y = "Posterior Membership Probability") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  theme(axis.title.y = element_text(face = "bold", size = 16)) +
  theme(axis.text.y = element_text(face = "bold", size = 10)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(strip.text.x = element_text(face = "bold", size = 16)) +
  theme(legend.title = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(face = "bold", size = 14)) +
  theme(legend.position = "bottom")
prow_dp15_ld01_dapc_10_p

# saved as A4, landscape orientation







#### ++ FIND.CLUSTERS - BIC - loop through K's (1-10) ####

## first check to see how many PCs will need to be retained in the loop (how many to get to 95% of cumulative variance?). Just using diffNgroup so that it doesn't ask me to select the number of clusters. Although, in this case, the number is selects is 2 (even though BIC K = 1 is the lowest)
prow_dp15_ld01_gl_pcs<- find.clusters(prow_dp15_ld01_gl, choose.n.clust = FALSE, criterion = "diffNgroup")

## looks like 25 PCs. Select 25 PCs and 2 clusters (will not be using this variable for anything other than checking this)
## add 25 in the n.pca loop below

maxK <- 10
prow_dp15_ld01_mat <- matrix(nrow=10, ncol=maxK)
colnames(prow_dp15_ld01_mat) <- 1:ncol(prow_dp15_ld01_mat)
for(i in 1:nrow(prow_dp15_ld01_mat)) {
  prow_dp15_ld01_gl_grp<-find.clusters(prow_dp15_ld01_gl, n.pca = 25, choose.n.clust = FALSE, max.n.clust = maxK)
  prow_dp15_ld01_mat[i,] <- prow_dp15_ld01_gl_grp$Kstat
}

prow_dp15_ld01_mat_df <- melt(prow_dp15_ld01_mat)
colnames(prow_dp15_ld01_mat_df)[1:3] <- c("Group", "K", "BIC")
prow_dp15_ld01_mat_df$K <- as.factor(prow_dp15_ld01_mat_df$K)
head(prow_dp15_ld01_mat_df)

ggplot(prow_dp15_ld01_mat_df, aes(x=K, y=BIC)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("A. protonotaria DP15, ld01") +
  xlab("Number of groups (K)")

# saved as A4, landscape orientation


# The best K = 1, but have to choose 2 for DAPC analysis or else it will fail.

# Do the same thing again but choose K = 2 (n.clust) and check the assignment plots. Keep n.pca as the same as the previous.

prow_dp15_ld01_gl_grp<-find.clusters(prow_dp15_ld01_gl, n.pca = 25, n.clust = 2)

# See the group assignments of each sampling location (i.e.,'pop' assigned upon loading in data)
table(pop(prow_dp15_ld01_gl), prow_dp15_ld01_gl_grp$grp)

# Compare K = 2 assignment with pop

table.value(table(pop(prow_dp15_ld01_gl), prow_dp15_ld01_gl_grp$grp), col.lab=paste("inferred", 1:2),
            row.lab=paste("original", 1:5))




#### ++ FIND.CLUSTERS - DAPC - biallelic only, 30x cross-validation loop ####
set.seed(NULL)

prow_dp15_ld01_finddapc_cv_loopi_n.pca<-list()
prow_dp15_ld01_finddapc_cv_loopi_n.da<-list()

for (i in 1:100)
{
  print(i)
  prow_dp15_ld01_finddapc_cv_i <- xvalDapc(tab(prow_dp15_ld01_gl, NA.method = "mean"), grp = prow_dp15_ld01_gl_grp$grp, n.rep = 30, parallel = "multicore", ncpus = 6L)
  
  prow_dp15_ld01_finddapc_cv_loopi_n.pca[[i]] <- prow_dp15_ld01_finddapc_cv_i$DAPC$n.pca
  prow_dp15_ld01_finddapc_cv_loopi_n.da[[i]] <- prow_dp15_ld01_finddapc_cv_i$DAPC$n.da
}


# check some basic stats

max(unlist(prow_dp15_ld01_finddapc_cv_loopi_n.pca))
min(unlist(prow_dp15_ld01_finddapc_cv_loopi_n.pca))
mean(unlist(prow_dp15_ld01_finddapc_cv_loopi_n.pca))
median(unlist(prow_dp15_ld01_finddapc_cv_loopi_n.pca))
hist(unlist(prow_dp15_ld01_finddapc_cv_loopi_n.pca))

# only run 'save' if needed to replace previous files
# save(prow_dp15_ld01_finddapc_cv_loopi_n.pca, file='./rda_files/prow_dp15_ld01_finddapc_cv_loopi_n.pca.rda')
# save(prow_dp15_ld01_finddapc_cv_loopi_n.da, file='./rda_files/prow_dp15_ld01_finddapc_cv_loopi_n.da.rda')



# reformat for ggplot histogram

prow_dp15_ld01_finddapc_cv_loopi_n.pca_df<-as.data.frame(matrix(unlist(prow_dp15_ld01_finddapc_cv_loopi_n.pca),nrow=length(prow_dp15_ld01_finddapc_cv_loopi_n.pca),byrow=TRUE))
colnames(prow_dp15_ld01_finddapc_cv_loopi_n.pca_df)[colnames(prow_dp15_ld01_finddapc_cv_loopi_n.pca_df) == 'V1'] <- 'Number_of_PCs'

# only run 'save' if needed to replace previous file
# save(prow_dp15_ld01_finddapc_cv_loopi_n.pca_df, file='./rda_files/prow_dp15_ld01_finddapc_cv_loopi_n.pca_df.rda')


ggplot(prow_dp15_ld01_finddapc_cv_loopi_n.pca_df, aes(x = Number_of_PCs)) +
  geom_bar() +
  scale_x_continuous(breaks=seq(2,25, by = 2)) +
  scale_y_continuous(breaks=seq(0,100, by = 5)) +
  ggtitle("Distribution of PCs after 100 cross-validations (30 reps) \nPROW DP15 LD01 using find.clusters") +
  xlab("Number of PCs retained") +
  ylab("Count")+
  theme_bw()

# Saved as A4 landscape 

# Looks like 1 is obviously our best bet, but let's try with more x-vals around the bulk of the distribution (1-2)






#### +++ 1000 xvals centered around # of PCs determined above ####

prow_dp15_ld01_finddapc_cv1000 <- xvalDapc(tab(prow_dp15_ld01_gl, NA.method = "mean"), grp = prow_dp15_ld01_gl_grp$grp, n.pca = 1:2, n.rep = 1000, parallel = "multicore", ncpus = 6L)

prow_dp15_ld01_finddapc_cv1000[-1] # 1 PCs

# $n.pca: 1 first PCs of PCA used
# $n.da: 1 discriminant functions saved
# $var (proportion of conserved variance): 0.051

# only run 'save' if needed to replace previous file
# save(prow_dp15_ld01_finddapc_cv1000, file='./rda_files/prow_dp15_ld01_finddapc_cv1000.rda')
load(file='./rda_files/prow_dp15_ld01_finddapc_cv1000.rda')


# This automatically puts the number of PCs determined in the xval
# legends
scatter.dapc(prow_dp15_ld01_finddapc_cv1000$DAPC, cex = 2, legend = TRUE,
             clabel = FALSE, cleg=1.5, posi.leg = "topright", scree.pca = TRUE,
             posi.pca = "topleft", posi.da="none", xax = 1, yax = 2, inset.solid = 1, col=c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), pch=19, bg.inset = "gray97")

# saved as A4, landscape orientation

# no legends
scatter.dapc(prow_dp15_ld01_finddapc_cv1000$DAPC, cex = 2, legend = FALSE,
             clabel = FALSE, cleg=1.5, scree.pca = FALSE,
             posi.da="none", xax = 1, yax = 2, inset.solid = 1, col=c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), pch=19, bg.inset = "gray97")

# saved as A4, landscape orientation

# assignment plot, shows the accuracy of the model to recover the each sample's sample location... cell color indicates the probability of assigning a given sample to the corresponding location (i.e., membership probability). Red = probability of 1 to white = probability of 0. Blue X's indicate original assigned location cluster.
assignplot(prow_dp15_ld01_finddapc_cv1000$DAPC, cex.lab=0.5, pch=4)

# saved as A4, landscape orientation


# here you can choose # of PCs to retain and use the n.da determined from above to make compoplot
prow_dp15_ld01_finddapc_1 <- dapc(prow_dp15_ld01_gl, n.pca = 1, n.da = 1)

compoplot(prow_dp15_ld01_finddapc_1, col = cols, posi = 'top')
prow_dp15_ld01_finddapc_1.results <- as.data.frame(prow_dp15_ld01_finddapc_1$posterior)
prow_dp15_ld01_finddapc_1.results$pop <- pop(prow_dp15_ld01_gl)
prow_dp15_ld01_finddapc_1.results$indNames <- rownames(prow_dp15_ld01_finddapc_1.results)
prow_dp15_ld01_finddapc_1.results <- pivot_longer(prow_dp15_ld01_finddapc_1.results, -c(pop, indNames))
colnames(prow_dp15_ld01_finddapc_1.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

prow_dp15_ld01_finddapc_1_p <- ggplot(prow_dp15_ld01_finddapc_1.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') +
  theme_few() +
  scale_fill_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3"), name = "Assigned\nPopulation") +
  facet_grid(~Original_Pop, scales = "free") +
  labs(y = "Posterior Membership Probability") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  theme(axis.title.y = element_text(face = "bold", size = 16)) +
  theme(axis.text.y = element_text(face = "bold", size = 10)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(strip.text.x = element_text(face = "bold", size = 16)) +
  theme(legend.title = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(face = "bold", size = 14)) +
  theme(legend.position = "bottom")
prow_dp15_ld01_finddapc_1_p

# saved as A4, landscape orientation

