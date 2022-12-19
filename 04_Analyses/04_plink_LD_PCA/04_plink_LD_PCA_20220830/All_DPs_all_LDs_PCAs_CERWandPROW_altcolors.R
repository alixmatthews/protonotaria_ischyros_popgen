#### LOAD LIBRARIES ####

setwd("~/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_snp/04_Analyses/04_plink_LD_PCA_20220830")

# following this tutorial: https://speciationgenomics.github.io/pca/

library(tidyverse)
library(ggplot2)
library(ggpubr)


#### CERW DEPTH SET AT 5X ####

#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld0001 ####
cerw_dp05_maxmiss100_ld0001_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld0001.eigenvec", col_names = FALSE)
cerw_dp05_maxmiss100_ld0001_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld0001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp05_maxmiss100_ld0001_pca <- cerw_dp05_maxmiss100_ld0001_pca[,-1]
# set names
names(cerw_dp05_maxmiss100_ld0001_pca)[1] <- "ind"

names(cerw_dp05_maxmiss100_ld0001_pca)[2:ncol(cerw_dp05_maxmiss100_ld0001_pca)] <- paste0("PC", 1:(ncol(cerw_dp05_maxmiss100_ld0001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp05_maxmiss100_ld0001_loc <- rep(NA, length(cerw_dp05_maxmiss100_ld0001_pca$ind))
cerw_dp05_maxmiss100_ld0001_loc[grep("IN", cerw_dp05_maxmiss100_ld0001_pca$ind)] <- "Indiana"
cerw_dp05_maxmiss100_ld0001_loc[grep("TN", cerw_dp05_maxmiss100_ld0001_pca$ind)] <- "Tennessee"
cerw_dp05_maxmiss100_ld0001_loc[grep("OZ", cerw_dp05_maxmiss100_ld0001_pca$ind)] <- "Ozarks"
cerw_dp05_maxmiss100_ld0001_loc[grep("WI", cerw_dp05_maxmiss100_ld0001_pca$ind)] <- "Wisconsin"
cerw_dp05_maxmiss100_ld0001_loc[grep("PA", cerw_dp05_maxmiss100_ld0001_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp05_maxmiss100_ld0001_pca <- as_tibble(data.frame(cerw_dp05_maxmiss100_ld0001_pca, cerw_dp05_maxmiss100_ld0001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp05_maxmiss100_ld0001_pve <- data.frame(PC = 1:20, pve = cerw_dp05_maxmiss100_ld0001_eigenval/sum(cerw_dp05_maxmiss100_ld0001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp05_maxmiss100_ld0001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp05_maxmiss100_ld0001_pve$pve)

# plot pca
cerw_dp05_maxmiss100_ld0001_pca_plot<-ggplot(cerw_dp05_maxmiss100_ld0001_pca, aes(PC1, PC2, col = cerw_dp05_maxmiss100_ld0001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.001)) +
    xlab(paste0("PC1 (", signif(cerw_dp05_maxmiss100_ld0001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp05_maxmiss100_ld0001_pve$pve[2], 3), "%)"))

cerw_dp05_maxmiss100_ld0001_pca_plot









#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld001 ####
cerw_dp05_maxmiss100_ld001_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld001.eigenvec", col_names = FALSE)
cerw_dp05_maxmiss100_ld001_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp05_maxmiss100_ld001_pca <- cerw_dp05_maxmiss100_ld001_pca[,-1]
# set names
names(cerw_dp05_maxmiss100_ld001_pca)[1] <- "ind"
names(cerw_dp05_maxmiss100_ld001_pca)[2:ncol(cerw_dp05_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(cerw_dp05_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp05_maxmiss100_ld001_loc <- rep(NA, length(cerw_dp05_maxmiss100_ld001_pca$ind))
cerw_dp05_maxmiss100_ld001_loc[grep("IN", cerw_dp05_maxmiss100_ld001_pca$ind)] <- "Indiana"
cerw_dp05_maxmiss100_ld001_loc[grep("TN", cerw_dp05_maxmiss100_ld001_pca$ind)] <- "Tennessee"
cerw_dp05_maxmiss100_ld001_loc[grep("OZ", cerw_dp05_maxmiss100_ld001_pca$ind)] <- "Ozarks"
cerw_dp05_maxmiss100_ld001_loc[grep("WI", cerw_dp05_maxmiss100_ld001_pca$ind)] <- "Wisconsin"
cerw_dp05_maxmiss100_ld001_loc[grep("PA", cerw_dp05_maxmiss100_ld001_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp05_maxmiss100_ld001_pca <- as_tibble(data.frame(cerw_dp05_maxmiss100_ld001_pca, cerw_dp05_maxmiss100_ld001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp05_maxmiss100_ld001_pve <- data.frame(PC = 1:20, pve = cerw_dp05_maxmiss100_ld001_eigenval/sum(cerw_dp05_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp05_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp05_maxmiss100_ld001_pve$pve)

# plot pca
cerw_dp05_maxmiss100_ld001_pca_plot<-ggplot(cerw_dp05_maxmiss100_ld001_pca, aes(PC1, PC2, col = cerw_dp05_maxmiss100_ld001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(cerw_dp05_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp05_maxmiss100_ld001_pve$pve[2], 3), "%)"))

cerw_dp05_maxmiss100_ld001_pca_plot















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld01 ####
cerw_dp05_maxmiss100_ld01_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld01.eigenvec", col_names = FALSE)
cerw_dp05_maxmiss100_ld01_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld01.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp05_maxmiss100_ld01_pca <- cerw_dp05_maxmiss100_ld01_pca[,-1]
# set names
names(cerw_dp05_maxmiss100_ld01_pca)[1] <- "ind"
names(cerw_dp05_maxmiss100_ld01_pca)[2:ncol(cerw_dp05_maxmiss100_ld01_pca)] <- paste0("PC", 1:(ncol(cerw_dp05_maxmiss100_ld01_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp05_maxmiss100_ld01_loc <- rep(NA, length(cerw_dp05_maxmiss100_ld01_pca$ind))
cerw_dp05_maxmiss100_ld01_loc[grep("IN", cerw_dp05_maxmiss100_ld01_pca$ind)] <- "Indiana"
cerw_dp05_maxmiss100_ld01_loc[grep("TN", cerw_dp05_maxmiss100_ld01_pca$ind)] <- "Tennessee"
cerw_dp05_maxmiss100_ld01_loc[grep("OZ", cerw_dp05_maxmiss100_ld01_pca$ind)] <- "Ozarks"
cerw_dp05_maxmiss100_ld01_loc[grep("WI", cerw_dp05_maxmiss100_ld01_pca$ind)] <- "Wisconsin"
cerw_dp05_maxmiss100_ld01_loc[grep("PA", cerw_dp05_maxmiss100_ld01_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp05_maxmiss100_ld01_pca <- as_tibble(data.frame(cerw_dp05_maxmiss100_ld01_pca, cerw_dp05_maxmiss100_ld01_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp05_maxmiss100_ld01_pve <- data.frame(PC = 1:20, pve = cerw_dp05_maxmiss100_ld01_eigenval/sum(cerw_dp05_maxmiss100_ld01_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp05_maxmiss100_ld01_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp05_maxmiss100_ld01_pve$pve)

# plot pca
cerw_dp05_maxmiss100_ld01_pca_plot<-ggplot(cerw_dp05_maxmiss100_ld01_pca, aes(PC1, PC2, col = cerw_dp05_maxmiss100_ld01_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.1)) +
    xlab(paste0("PC1 (", signif(cerw_dp05_maxmiss100_ld01_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp05_maxmiss100_ld01_pve$pve[2], 3), "%)"))

cerw_dp05_maxmiss100_ld01_pca_plot



#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld02 ####
cerw_dp05_maxmiss100_ld02_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld02.eigenvec", col_names = FALSE)
cerw_dp05_maxmiss100_ld02_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld02.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp05_maxmiss100_ld02_pca <- cerw_dp05_maxmiss100_ld02_pca[,-1]
# set names
names(cerw_dp05_maxmiss100_ld02_pca)[1] <- "ind"
names(cerw_dp05_maxmiss100_ld02_pca)[2:ncol(cerw_dp05_maxmiss100_ld02_pca)] <- paste0("PC", 1:(ncol(cerw_dp05_maxmiss100_ld02_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp05_maxmiss100_ld02_loc <- rep(NA, length(cerw_dp05_maxmiss100_ld02_pca$ind))
cerw_dp05_maxmiss100_ld02_loc[grep("IN", cerw_dp05_maxmiss100_ld02_pca$ind)] <- "Indiana"
cerw_dp05_maxmiss100_ld02_loc[grep("TN", cerw_dp05_maxmiss100_ld02_pca$ind)] <- "Tennessee"
cerw_dp05_maxmiss100_ld02_loc[grep("OZ", cerw_dp05_maxmiss100_ld02_pca$ind)] <- "Ozarks"
cerw_dp05_maxmiss100_ld02_loc[grep("WI", cerw_dp05_maxmiss100_ld02_pca$ind)] <- "Wisconsin"
cerw_dp05_maxmiss100_ld02_loc[grep("PA", cerw_dp05_maxmiss100_ld02_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp05_maxmiss100_ld02_pca <- as_tibble(data.frame(cerw_dp05_maxmiss100_ld02_pca, cerw_dp05_maxmiss100_ld02_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp05_maxmiss100_ld02_pve <- data.frame(PC = 1:20, pve = cerw_dp05_maxmiss100_ld02_eigenval/sum(cerw_dp05_maxmiss100_ld02_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp05_maxmiss100_ld02_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp05_maxmiss100_ld02_pve$pve)

# plot pca
cerw_dp05_maxmiss100_ld02_pca_plot<-ggplot(cerw_dp05_maxmiss100_ld02_pca, aes(PC1, PC2, col = cerw_dp05_maxmiss100_ld02_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.2)) +
    xlab(paste0("PC1 (", signif(cerw_dp05_maxmiss100_ld02_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp05_maxmiss100_ld02_pve$pve[2], 3), "%)"))

cerw_dp05_maxmiss100_ld02_pca_plot


















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld03 ####
cerw_dp05_maxmiss100_ld03_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld03.eigenvec", col_names = FALSE)
cerw_dp05_maxmiss100_ld03_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld03.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp05_maxmiss100_ld03_pca <- cerw_dp05_maxmiss100_ld03_pca[,-1]
# set names
names(cerw_dp05_maxmiss100_ld03_pca)[1] <- "ind"
names(cerw_dp05_maxmiss100_ld03_pca)[2:ncol(cerw_dp05_maxmiss100_ld03_pca)] <- paste0("PC", 1:(ncol(cerw_dp05_maxmiss100_ld03_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp05_maxmiss100_ld03_loc <- rep(NA, length(cerw_dp05_maxmiss100_ld03_pca$ind))
cerw_dp05_maxmiss100_ld03_loc[grep("IN", cerw_dp05_maxmiss100_ld03_pca$ind)] <- "Indiana"
cerw_dp05_maxmiss100_ld03_loc[grep("TN", cerw_dp05_maxmiss100_ld03_pca$ind)] <- "Tennessee"
cerw_dp05_maxmiss100_ld03_loc[grep("OZ", cerw_dp05_maxmiss100_ld03_pca$ind)] <- "Ozarks"
cerw_dp05_maxmiss100_ld03_loc[grep("WI", cerw_dp05_maxmiss100_ld03_pca$ind)] <- "Wisconsin"
cerw_dp05_maxmiss100_ld03_loc[grep("PA", cerw_dp05_maxmiss100_ld03_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp05_maxmiss100_ld03_pca <- as_tibble(data.frame(cerw_dp05_maxmiss100_ld03_pca, cerw_dp05_maxmiss100_ld03_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp05_maxmiss100_ld03_pve <- data.frame(PC = 1:20, pve = cerw_dp05_maxmiss100_ld03_eigenval/sum(cerw_dp05_maxmiss100_ld03_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp05_maxmiss100_ld03_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp05_maxmiss100_ld03_pve$pve)

# plot pca
cerw_dp05_maxmiss100_ld03_pca_plot<-ggplot(cerw_dp05_maxmiss100_ld03_pca, aes(PC1, PC2, col = cerw_dp05_maxmiss100_ld03_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.3)) +
    xlab(paste0("PC1 (", signif(cerw_dp05_maxmiss100_ld03_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp05_maxmiss100_ld03_pve$pve[2], 3), "%)"))

cerw_dp05_maxmiss100_ld03_pca_plot















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld04 ####
cerw_dp05_maxmiss100_ld04_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld04.eigenvec", col_names = FALSE)
cerw_dp05_maxmiss100_ld04_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld04.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp05_maxmiss100_ld04_pca <- cerw_dp05_maxmiss100_ld04_pca[,-1]
# set names
names(cerw_dp05_maxmiss100_ld04_pca)[1] <- "ind"
names(cerw_dp05_maxmiss100_ld04_pca)[2:ncol(cerw_dp05_maxmiss100_ld04_pca)] <- paste0("PC", 1:(ncol(cerw_dp05_maxmiss100_ld04_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp05_maxmiss100_ld04_loc <- rep(NA, length(cerw_dp05_maxmiss100_ld04_pca$ind))
cerw_dp05_maxmiss100_ld04_loc[grep("IN", cerw_dp05_maxmiss100_ld04_pca$ind)] <- "Indiana"
cerw_dp05_maxmiss100_ld04_loc[grep("TN", cerw_dp05_maxmiss100_ld04_pca$ind)] <- "Tennessee"
cerw_dp05_maxmiss100_ld04_loc[grep("OZ", cerw_dp05_maxmiss100_ld04_pca$ind)] <- "Ozarks"
cerw_dp05_maxmiss100_ld04_loc[grep("WI", cerw_dp05_maxmiss100_ld04_pca$ind)] <- "Wisconsin"
cerw_dp05_maxmiss100_ld04_loc[grep("PA", cerw_dp05_maxmiss100_ld04_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp05_maxmiss100_ld04_pca <- as_tibble(data.frame(cerw_dp05_maxmiss100_ld04_pca, cerw_dp05_maxmiss100_ld04_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp05_maxmiss100_ld04_pve <- data.frame(PC = 1:20, pve = cerw_dp05_maxmiss100_ld04_eigenval/sum(cerw_dp05_maxmiss100_ld04_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp05_maxmiss100_ld04_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp05_maxmiss100_ld04_pve$pve)

# plot pca
cerw_dp05_maxmiss100_ld04_pca_plot<-ggplot(cerw_dp05_maxmiss100_ld04_pca, aes(PC1, PC2, col = cerw_dp05_maxmiss100_ld04_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.4)) +
    xlab(paste0("PC1 (", signif(cerw_dp05_maxmiss100_ld04_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp05_maxmiss100_ld04_pve$pve[2], 3), "%)"))

cerw_dp05_maxmiss100_ld04_pca_plot




#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld05 ####
cerw_dp05_maxmiss100_ld05_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld05.eigenvec", col_names = FALSE)
cerw_dp05_maxmiss100_ld05_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp05_maxmiss100_ld05_pca <- cerw_dp05_maxmiss100_ld05_pca[,-1]
# set names
names(cerw_dp05_maxmiss100_ld05_pca)[1] <- "ind"
names(cerw_dp05_maxmiss100_ld05_pca)[2:ncol(cerw_dp05_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(cerw_dp05_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp05_maxmiss100_ld05_loc <- rep(NA, length(cerw_dp05_maxmiss100_ld05_pca$ind))
cerw_dp05_maxmiss100_ld05_loc[grep("IN", cerw_dp05_maxmiss100_ld05_pca$ind)] <- "Indiana"
cerw_dp05_maxmiss100_ld05_loc[grep("TN", cerw_dp05_maxmiss100_ld05_pca$ind)] <- "Tennessee"
cerw_dp05_maxmiss100_ld05_loc[grep("OZ", cerw_dp05_maxmiss100_ld05_pca$ind)] <- "Ozarks"
cerw_dp05_maxmiss100_ld05_loc[grep("WI", cerw_dp05_maxmiss100_ld05_pca$ind)] <- "Wisconsin"
cerw_dp05_maxmiss100_ld05_loc[grep("PA", cerw_dp05_maxmiss100_ld05_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp05_maxmiss100_ld05_pca <- as_tibble(data.frame(cerw_dp05_maxmiss100_ld05_pca, cerw_dp05_maxmiss100_ld05_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp05_maxmiss100_ld05_pve <- data.frame(PC = 1:20, pve = cerw_dp05_maxmiss100_ld05_eigenval/sum(cerw_dp05_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp05_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp05_maxmiss100_ld05_pve$pve)

# plot pca
cerw_dp05_maxmiss100_ld05_pca_plot<-ggplot(cerw_dp05_maxmiss100_ld05_pca, aes(PC1, PC2, col = cerw_dp05_maxmiss100_ld05_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(cerw_dp05_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp05_maxmiss100_ld05_pve$pve[2], 3), "%)"))

cerw_dp05_maxmiss100_ld05_pca_plot








#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld08 ####
cerw_dp05_maxmiss100_ld08_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld08.eigenvec", col_names = FALSE)
cerw_dp05_maxmiss100_ld08_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld08.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp05_maxmiss100_ld08_pca <- cerw_dp05_maxmiss100_ld08_pca[,-1]
# set names
names(cerw_dp05_maxmiss100_ld08_pca)[1] <- "ind"
names(cerw_dp05_maxmiss100_ld08_pca)[2:ncol(cerw_dp05_maxmiss100_ld08_pca)] <- paste0("PC", 1:(ncol(cerw_dp05_maxmiss100_ld08_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp05_maxmiss100_ld08_loc <- rep(NA, length(cerw_dp05_maxmiss100_ld08_pca$ind))
cerw_dp05_maxmiss100_ld08_loc[grep("IN", cerw_dp05_maxmiss100_ld08_pca$ind)] <- "Indiana"
cerw_dp05_maxmiss100_ld08_loc[grep("TN", cerw_dp05_maxmiss100_ld08_pca$ind)] <- "Tennessee"
cerw_dp05_maxmiss100_ld08_loc[grep("OZ", cerw_dp05_maxmiss100_ld08_pca$ind)] <- "Ozarks"
cerw_dp05_maxmiss100_ld08_loc[grep("WI", cerw_dp05_maxmiss100_ld08_pca$ind)] <- "Wisconsin"
cerw_dp05_maxmiss100_ld08_loc[grep("PA", cerw_dp05_maxmiss100_ld08_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp05_maxmiss100_ld08_pca <- as_tibble(data.frame(cerw_dp05_maxmiss100_ld08_pca, cerw_dp05_maxmiss100_ld08_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp05_maxmiss100_ld08_pve <- data.frame(PC = 1:20, pve = cerw_dp05_maxmiss100_ld08_eigenval/sum(cerw_dp05_maxmiss100_ld08_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp05_maxmiss100_ld08_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp05_maxmiss100_ld08_pve$pve)

# plot pca
cerw_dp05_maxmiss100_ld08_pca_plot<-ggplot(cerw_dp05_maxmiss100_ld08_pca, aes(PC1, PC2, col = cerw_dp05_maxmiss100_ld08_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.8)) +
    xlab(paste0("PC1 (", signif(cerw_dp05_maxmiss100_ld08_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp05_maxmiss100_ld08_pve$pve[2], 3), "%)"))

cerw_dp05_maxmiss100_ld08_pca_plot



















#### CERW DEPTH (DP) SET AT 10X ####


#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld0001 ####
cerw_dp10_maxmiss100_ld0001_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld0001.eigenvec", col_names = FALSE)
cerw_dp10_maxmiss100_ld0001_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld0001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp10_maxmiss100_ld0001_pca <- cerw_dp10_maxmiss100_ld0001_pca[,-1]
# set names
names(cerw_dp10_maxmiss100_ld0001_pca)[1] <- "ind"

names(cerw_dp10_maxmiss100_ld0001_pca)[2:ncol(cerw_dp10_maxmiss100_ld0001_pca)] <- paste0("PC", 1:(ncol(cerw_dp10_maxmiss100_ld0001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp10_maxmiss100_ld0001_loc <- rep(NA, length(cerw_dp10_maxmiss100_ld0001_pca$ind))
cerw_dp10_maxmiss100_ld0001_loc[grep("IN", cerw_dp10_maxmiss100_ld0001_pca$ind)] <- "Indiana"
cerw_dp10_maxmiss100_ld0001_loc[grep("TN", cerw_dp10_maxmiss100_ld0001_pca$ind)] <- "Tennessee"
cerw_dp10_maxmiss100_ld0001_loc[grep("OZ", cerw_dp10_maxmiss100_ld0001_pca$ind)] <- "Ozarks"
cerw_dp10_maxmiss100_ld0001_loc[grep("WI", cerw_dp10_maxmiss100_ld0001_pca$ind)] <- "Wisconsin"
cerw_dp10_maxmiss100_ld0001_loc[grep("PA", cerw_dp10_maxmiss100_ld0001_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp10_maxmiss100_ld0001_pca <- as_tibble(data.frame(cerw_dp10_maxmiss100_ld0001_pca, cerw_dp10_maxmiss100_ld0001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp10_maxmiss100_ld0001_pve <- data.frame(PC = 1:20, pve = cerw_dp10_maxmiss100_ld0001_eigenval/sum(cerw_dp10_maxmiss100_ld0001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp10_maxmiss100_ld0001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp10_maxmiss100_ld0001_pve$pve)

# plot pca
cerw_dp10_maxmiss100_ld0001_pca_plot<-ggplot(cerw_dp10_maxmiss100_ld0001_pca, aes(PC1, PC2, col = cerw_dp10_maxmiss100_ld0001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.001)) +
    xlab(paste0("PC1 (", signif(cerw_dp10_maxmiss100_ld0001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp10_maxmiss100_ld0001_pve$pve[2], 3), "%)"))

cerw_dp10_maxmiss100_ld0001_pca_plot











#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001 ####
cerw_dp10_maxmiss100_ld001_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.eigenvec", col_names = FALSE)
cerw_dp10_maxmiss100_ld001_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp10_maxmiss100_ld001_pca <- cerw_dp10_maxmiss100_ld001_pca[,-1]
# set names
names(cerw_dp10_maxmiss100_ld001_pca)[1] <- "ind"
names(cerw_dp10_maxmiss100_ld001_pca)[2:ncol(cerw_dp10_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(cerw_dp10_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp10_maxmiss100_ld001_loc <- rep(NA, length(cerw_dp10_maxmiss100_ld001_pca$ind))
cerw_dp10_maxmiss100_ld001_loc[grep("IN", cerw_dp10_maxmiss100_ld001_pca$ind)] <- "Indiana"
cerw_dp10_maxmiss100_ld001_loc[grep("TN", cerw_dp10_maxmiss100_ld001_pca$ind)] <- "Tennessee"
cerw_dp10_maxmiss100_ld001_loc[grep("OZ", cerw_dp10_maxmiss100_ld001_pca$ind)] <- "Ozarks"
cerw_dp10_maxmiss100_ld001_loc[grep("WI", cerw_dp10_maxmiss100_ld001_pca$ind)] <- "Wisconsin"
cerw_dp10_maxmiss100_ld001_loc[grep("PA", cerw_dp10_maxmiss100_ld001_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp10_maxmiss100_ld001_pca <- as_tibble(data.frame(cerw_dp10_maxmiss100_ld001_pca, cerw_dp10_maxmiss100_ld001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp10_maxmiss100_ld001_pve <- data.frame(PC = 1:20, pve = cerw_dp10_maxmiss100_ld001_eigenval/sum(cerw_dp10_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp10_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp10_maxmiss100_ld001_pve$pve)

# plot pca
cerw_dp10_maxmiss100_ld001_pca_plot<-ggplot(cerw_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = cerw_dp10_maxmiss100_ld001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(cerw_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)"))

cerw_dp10_maxmiss100_ld001_pca_plot















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld01 ####
cerw_dp10_maxmiss100_ld01_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld01.eigenvec", col_names = FALSE)
cerw_dp10_maxmiss100_ld01_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld01.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp10_maxmiss100_ld01_pca <- cerw_dp10_maxmiss100_ld01_pca[,-1]
# set names
names(cerw_dp10_maxmiss100_ld01_pca)[1] <- "ind"
names(cerw_dp10_maxmiss100_ld01_pca)[2:ncol(cerw_dp10_maxmiss100_ld01_pca)] <- paste0("PC", 1:(ncol(cerw_dp10_maxmiss100_ld01_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp10_maxmiss100_ld01_loc <- rep(NA, length(cerw_dp10_maxmiss100_ld01_pca$ind))
cerw_dp10_maxmiss100_ld01_loc[grep("IN", cerw_dp10_maxmiss100_ld01_pca$ind)] <- "Indiana"
cerw_dp10_maxmiss100_ld01_loc[grep("TN", cerw_dp10_maxmiss100_ld01_pca$ind)] <- "Tennessee"
cerw_dp10_maxmiss100_ld01_loc[grep("OZ", cerw_dp10_maxmiss100_ld01_pca$ind)] <- "Ozarks"
cerw_dp10_maxmiss100_ld01_loc[grep("WI", cerw_dp10_maxmiss100_ld01_pca$ind)] <- "Wisconsin"
cerw_dp10_maxmiss100_ld01_loc[grep("PA", cerw_dp10_maxmiss100_ld01_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp10_maxmiss100_ld01_pca <- as_tibble(data.frame(cerw_dp10_maxmiss100_ld01_pca, cerw_dp10_maxmiss100_ld01_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp10_maxmiss100_ld01_pve <- data.frame(PC = 1:20, pve = cerw_dp10_maxmiss100_ld01_eigenval/sum(cerw_dp10_maxmiss100_ld01_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp10_maxmiss100_ld01_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp10_maxmiss100_ld01_pve$pve)

# plot pca
cerw_dp10_maxmiss100_ld01_pca_plot<-ggplot(cerw_dp10_maxmiss100_ld01_pca, aes(PC1, PC2, col = cerw_dp10_maxmiss100_ld01_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.1)) +
    xlab(paste0("PC1 (", signif(cerw_dp10_maxmiss100_ld01_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp10_maxmiss100_ld01_pve$pve[2], 3), "%)"))

cerw_dp10_maxmiss100_ld01_pca_plot



#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld02 ####
cerw_dp10_maxmiss100_ld02_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld02.eigenvec", col_names = FALSE)
cerw_dp10_maxmiss100_ld02_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld02.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp10_maxmiss100_ld02_pca <- cerw_dp10_maxmiss100_ld02_pca[,-1]
# set names
names(cerw_dp10_maxmiss100_ld02_pca)[1] <- "ind"
names(cerw_dp10_maxmiss100_ld02_pca)[2:ncol(cerw_dp10_maxmiss100_ld02_pca)] <- paste0("PC", 1:(ncol(cerw_dp10_maxmiss100_ld02_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp10_maxmiss100_ld02_loc <- rep(NA, length(cerw_dp10_maxmiss100_ld02_pca$ind))
cerw_dp10_maxmiss100_ld02_loc[grep("IN", cerw_dp10_maxmiss100_ld02_pca$ind)] <- "Indiana"
cerw_dp10_maxmiss100_ld02_loc[grep("TN", cerw_dp10_maxmiss100_ld02_pca$ind)] <- "Tennessee"
cerw_dp10_maxmiss100_ld02_loc[grep("OZ", cerw_dp10_maxmiss100_ld02_pca$ind)] <- "Ozarks"
cerw_dp10_maxmiss100_ld02_loc[grep("WI", cerw_dp10_maxmiss100_ld02_pca$ind)] <- "Wisconsin"
cerw_dp10_maxmiss100_ld02_loc[grep("PA", cerw_dp10_maxmiss100_ld02_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp10_maxmiss100_ld02_pca <- as_tibble(data.frame(cerw_dp10_maxmiss100_ld02_pca, cerw_dp10_maxmiss100_ld02_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp10_maxmiss100_ld02_pve <- data.frame(PC = 1:20, pve = cerw_dp10_maxmiss100_ld02_eigenval/sum(cerw_dp10_maxmiss100_ld02_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp10_maxmiss100_ld02_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp10_maxmiss100_ld02_pve$pve)

# plot pca
cerw_dp10_maxmiss100_ld02_pca_plot<-ggplot(cerw_dp10_maxmiss100_ld02_pca, aes(PC1, PC2, col = cerw_dp10_maxmiss100_ld02_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.2)) +
    xlab(paste0("PC1 (", signif(cerw_dp10_maxmiss100_ld02_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp10_maxmiss100_ld02_pve$pve[2], 3), "%)"))

cerw_dp10_maxmiss100_ld02_pca_plot


















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld03 ####
cerw_dp10_maxmiss100_ld03_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld03.eigenvec", col_names = FALSE)
cerw_dp10_maxmiss100_ld03_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld03.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp10_maxmiss100_ld03_pca <- cerw_dp10_maxmiss100_ld03_pca[,-1]
# set names
names(cerw_dp10_maxmiss100_ld03_pca)[1] <- "ind"
names(cerw_dp10_maxmiss100_ld03_pca)[2:ncol(cerw_dp10_maxmiss100_ld03_pca)] <- paste0("PC", 1:(ncol(cerw_dp10_maxmiss100_ld03_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp10_maxmiss100_ld03_loc <- rep(NA, length(cerw_dp10_maxmiss100_ld03_pca$ind))
cerw_dp10_maxmiss100_ld03_loc[grep("IN", cerw_dp10_maxmiss100_ld03_pca$ind)] <- "Indiana"
cerw_dp10_maxmiss100_ld03_loc[grep("TN", cerw_dp10_maxmiss100_ld03_pca$ind)] <- "Tennessee"
cerw_dp10_maxmiss100_ld03_loc[grep("OZ", cerw_dp10_maxmiss100_ld03_pca$ind)] <- "Ozarks"
cerw_dp10_maxmiss100_ld03_loc[grep("WI", cerw_dp10_maxmiss100_ld03_pca$ind)] <- "Wisconsin"
cerw_dp10_maxmiss100_ld03_loc[grep("PA", cerw_dp10_maxmiss100_ld03_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp10_maxmiss100_ld03_pca <- as_tibble(data.frame(cerw_dp10_maxmiss100_ld03_pca, cerw_dp10_maxmiss100_ld03_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp10_maxmiss100_ld03_pve <- data.frame(PC = 1:20, pve = cerw_dp10_maxmiss100_ld03_eigenval/sum(cerw_dp10_maxmiss100_ld03_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp10_maxmiss100_ld03_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp10_maxmiss100_ld03_pve$pve)

# plot pca
cerw_dp10_maxmiss100_ld03_pca_plot<-ggplot(cerw_dp10_maxmiss100_ld03_pca, aes(PC1, PC2, col = cerw_dp10_maxmiss100_ld03_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.3)) +
    xlab(paste0("PC1 (", signif(cerw_dp10_maxmiss100_ld03_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp10_maxmiss100_ld03_pve$pve[2], 3), "%)"))

cerw_dp10_maxmiss100_ld03_pca_plot















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld04 ####
cerw_dp10_maxmiss100_ld04_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld04.eigenvec", col_names = FALSE)
cerw_dp10_maxmiss100_ld04_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld04.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp10_maxmiss100_ld04_pca <- cerw_dp10_maxmiss100_ld04_pca[,-1]
# set names
names(cerw_dp10_maxmiss100_ld04_pca)[1] <- "ind"
names(cerw_dp10_maxmiss100_ld04_pca)[2:ncol(cerw_dp10_maxmiss100_ld04_pca)] <- paste0("PC", 1:(ncol(cerw_dp10_maxmiss100_ld04_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp10_maxmiss100_ld04_loc <- rep(NA, length(cerw_dp10_maxmiss100_ld04_pca$ind))
cerw_dp10_maxmiss100_ld04_loc[grep("IN", cerw_dp10_maxmiss100_ld04_pca$ind)] <- "Indiana"
cerw_dp10_maxmiss100_ld04_loc[grep("TN", cerw_dp10_maxmiss100_ld04_pca$ind)] <- "Tennessee"
cerw_dp10_maxmiss100_ld04_loc[grep("OZ", cerw_dp10_maxmiss100_ld04_pca$ind)] <- "Ozarks"
cerw_dp10_maxmiss100_ld04_loc[grep("WI", cerw_dp10_maxmiss100_ld04_pca$ind)] <- "Wisconsin"
cerw_dp10_maxmiss100_ld04_loc[grep("PA", cerw_dp10_maxmiss100_ld04_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp10_maxmiss100_ld04_pca <- as_tibble(data.frame(cerw_dp10_maxmiss100_ld04_pca, cerw_dp10_maxmiss100_ld04_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp10_maxmiss100_ld04_pve <- data.frame(PC = 1:20, pve = cerw_dp10_maxmiss100_ld04_eigenval/sum(cerw_dp10_maxmiss100_ld04_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp10_maxmiss100_ld04_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp10_maxmiss100_ld04_pve$pve)

# plot pca
cerw_dp10_maxmiss100_ld04_pca_plot<-ggplot(cerw_dp10_maxmiss100_ld04_pca, aes(PC1, PC2, col = cerw_dp10_maxmiss100_ld04_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.4)) +
    xlab(paste0("PC1 (", signif(cerw_dp10_maxmiss100_ld04_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp10_maxmiss100_ld04_pve$pve[2], 3), "%)"))

cerw_dp10_maxmiss100_ld04_pca_plot




#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld05 ####
cerw_dp10_maxmiss100_ld05_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld05.eigenvec", col_names = FALSE)
cerw_dp10_maxmiss100_ld05_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp10_maxmiss100_ld05_pca <- cerw_dp10_maxmiss100_ld05_pca[,-1]
# set names
names(cerw_dp10_maxmiss100_ld05_pca)[1] <- "ind"
names(cerw_dp10_maxmiss100_ld05_pca)[2:ncol(cerw_dp10_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(cerw_dp10_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp10_maxmiss100_ld05_loc <- rep(NA, length(cerw_dp10_maxmiss100_ld05_pca$ind))
cerw_dp10_maxmiss100_ld05_loc[grep("IN", cerw_dp10_maxmiss100_ld05_pca$ind)] <- "Indiana"
cerw_dp10_maxmiss100_ld05_loc[grep("TN", cerw_dp10_maxmiss100_ld05_pca$ind)] <- "Tennessee"
cerw_dp10_maxmiss100_ld05_loc[grep("OZ", cerw_dp10_maxmiss100_ld05_pca$ind)] <- "Ozarks"
cerw_dp10_maxmiss100_ld05_loc[grep("WI", cerw_dp10_maxmiss100_ld05_pca$ind)] <- "Wisconsin"
cerw_dp10_maxmiss100_ld05_loc[grep("PA", cerw_dp10_maxmiss100_ld05_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp10_maxmiss100_ld05_pca <- as_tibble(data.frame(cerw_dp10_maxmiss100_ld05_pca, cerw_dp10_maxmiss100_ld05_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp10_maxmiss100_ld05_pve <- data.frame(PC = 1:20, pve = cerw_dp10_maxmiss100_ld05_eigenval/sum(cerw_dp10_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp10_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp10_maxmiss100_ld05_pve$pve)

# plot pca
cerw_dp10_maxmiss100_ld05_pca_plot<-ggplot(cerw_dp10_maxmiss100_ld05_pca, aes(PC1, PC2, col = cerw_dp10_maxmiss100_ld05_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(cerw_dp10_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp10_maxmiss100_ld05_pve$pve[2], 3), "%)"))

cerw_dp10_maxmiss100_ld05_pca_plot














#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld08 ####
cerw_dp10_maxmiss100_ld08_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld08.eigenvec", col_names = FALSE)
cerw_dp10_maxmiss100_ld08_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld08.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp10_maxmiss100_ld08_pca <- cerw_dp10_maxmiss100_ld08_pca[,-1]
# set names
names(cerw_dp10_maxmiss100_ld08_pca)[1] <- "ind"
names(cerw_dp10_maxmiss100_ld08_pca)[2:ncol(cerw_dp10_maxmiss100_ld08_pca)] <- paste0("PC", 1:(ncol(cerw_dp10_maxmiss100_ld08_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp10_maxmiss100_ld08_loc <- rep(NA, length(cerw_dp10_maxmiss100_ld08_pca$ind))
cerw_dp10_maxmiss100_ld08_loc[grep("IN", cerw_dp10_maxmiss100_ld08_pca$ind)] <- "Indiana"
cerw_dp10_maxmiss100_ld08_loc[grep("TN", cerw_dp10_maxmiss100_ld08_pca$ind)] <- "Tennessee"
cerw_dp10_maxmiss100_ld08_loc[grep("OZ", cerw_dp10_maxmiss100_ld08_pca$ind)] <- "Ozarks"
cerw_dp10_maxmiss100_ld08_loc[grep("WI", cerw_dp10_maxmiss100_ld08_pca$ind)] <- "Wisconsin"
cerw_dp10_maxmiss100_ld08_loc[grep("PA", cerw_dp10_maxmiss100_ld08_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp10_maxmiss100_ld08_pca <- as_tibble(data.frame(cerw_dp10_maxmiss100_ld08_pca, cerw_dp10_maxmiss100_ld08_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp10_maxmiss100_ld08_pve <- data.frame(PC = 1:20, pve = cerw_dp10_maxmiss100_ld08_eigenval/sum(cerw_dp10_maxmiss100_ld08_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp10_maxmiss100_ld08_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp10_maxmiss100_ld08_pve$pve)

# plot pca
cerw_dp10_maxmiss100_ld08_pca_plot<-ggplot(cerw_dp10_maxmiss100_ld08_pca, aes(PC1, PC2, col = cerw_dp10_maxmiss100_ld08_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.8)) +
    xlab(paste0("PC1 (", signif(cerw_dp10_maxmiss100_ld08_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp10_maxmiss100_ld08_pve$pve[2], 3), "%)"))

cerw_dp10_maxmiss100_ld08_pca_plot













#### CERW DEPTH (DP) SET AT 15X ####

#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld0001 ####
cerw_dp15_maxmiss100_ld0001_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld0001.eigenvec", col_names = FALSE)
cerw_dp15_maxmiss100_ld0001_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld0001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp15_maxmiss100_ld0001_pca <- cerw_dp15_maxmiss100_ld0001_pca[,-1]
# set names
names(cerw_dp15_maxmiss100_ld0001_pca)[1] <- "ind"
names(cerw_dp15_maxmiss100_ld0001_pca)[2:ncol(cerw_dp15_maxmiss100_ld0001_pca)] <- paste0("PC", 1:(ncol(cerw_dp15_maxmiss100_ld0001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp15_maxmiss100_ld0001_loc <- rep(NA, length(cerw_dp15_maxmiss100_ld0001_pca$ind))
cerw_dp15_maxmiss100_ld0001_loc[grep("IN", cerw_dp15_maxmiss100_ld0001_pca$ind)] <- "Indiana"
cerw_dp15_maxmiss100_ld0001_loc[grep("TN", cerw_dp15_maxmiss100_ld0001_pca$ind)] <- "Tennessee"
cerw_dp15_maxmiss100_ld0001_loc[grep("OZ", cerw_dp15_maxmiss100_ld0001_pca$ind)] <- "Ozarks"
cerw_dp15_maxmiss100_ld0001_loc[grep("WI", cerw_dp15_maxmiss100_ld0001_pca$ind)] <- "Wisconsin"
cerw_dp15_maxmiss100_ld0001_loc[grep("PA", cerw_dp15_maxmiss100_ld0001_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp15_maxmiss100_ld0001_pca <- as_tibble(data.frame(cerw_dp15_maxmiss100_ld0001_pca, cerw_dp15_maxmiss100_ld0001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp15_maxmiss100_ld0001_pve <- data.frame(PC = 1:20, pve = cerw_dp15_maxmiss100_ld0001_eigenval/sum(cerw_dp15_maxmiss100_ld0001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp15_maxmiss100_ld0001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp15_maxmiss100_ld0001_pve$pve)

# plot pca
cerw_dp15_maxmiss100_ld0001_pca_plot<-ggplot(cerw_dp15_maxmiss100_ld0001_pca, aes(PC1, PC2, col = cerw_dp15_maxmiss100_ld0001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.001)) +
    xlab(paste0("PC1 (", signif(cerw_dp15_maxmiss100_ld0001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp15_maxmiss100_ld0001_pve$pve[2], 3), "%)"))

cerw_dp15_maxmiss100_ld0001_pca_plot











#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001 ####
cerw_dp15_maxmiss100_ld001_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001.eigenvec", col_names = FALSE)
cerw_dp15_maxmiss100_ld001_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp15_maxmiss100_ld001_pca <- cerw_dp15_maxmiss100_ld001_pca[,-1]
# set names
names(cerw_dp15_maxmiss100_ld001_pca)[1] <- "ind"
names(cerw_dp15_maxmiss100_ld001_pca)[2:ncol(cerw_dp15_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(cerw_dp15_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp15_maxmiss100_ld001_loc <- rep(NA, length(cerw_dp15_maxmiss100_ld001_pca$ind))
cerw_dp15_maxmiss100_ld001_loc[grep("IN", cerw_dp15_maxmiss100_ld001_pca$ind)] <- "Indiana"
cerw_dp15_maxmiss100_ld001_loc[grep("TN", cerw_dp15_maxmiss100_ld001_pca$ind)] <- "Tennessee"
cerw_dp15_maxmiss100_ld001_loc[grep("OZ", cerw_dp15_maxmiss100_ld001_pca$ind)] <- "Ozarks"
cerw_dp15_maxmiss100_ld001_loc[grep("WI", cerw_dp15_maxmiss100_ld001_pca$ind)] <- "Wisconsin"
cerw_dp15_maxmiss100_ld001_loc[grep("PA", cerw_dp15_maxmiss100_ld001_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp15_maxmiss100_ld001_pca <- as_tibble(data.frame(cerw_dp15_maxmiss100_ld001_pca, cerw_dp15_maxmiss100_ld001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp15_maxmiss100_ld001_pve <- data.frame(PC = 1:20, pve = cerw_dp15_maxmiss100_ld001_eigenval/sum(cerw_dp15_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp15_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp15_maxmiss100_ld001_pve$pve)

# plot pca
cerw_dp15_maxmiss100_ld001_pca_plot<-ggplot(cerw_dp15_maxmiss100_ld001_pca, aes(PC1, PC2, col = cerw_dp15_maxmiss100_ld001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(cerw_dp15_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp15_maxmiss100_ld001_pve$pve[2], 3), "%)"))

cerw_dp15_maxmiss100_ld001_pca_plot















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01 ####
cerw_dp15_maxmiss100_ld01_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.eigenvec", col_names = FALSE)
cerw_dp15_maxmiss100_ld01_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp15_maxmiss100_ld01_pca <- cerw_dp15_maxmiss100_ld01_pca[,-1]
# set names
names(cerw_dp15_maxmiss100_ld01_pca)[1] <- "ind"
names(cerw_dp15_maxmiss100_ld01_pca)[2:ncol(cerw_dp15_maxmiss100_ld01_pca)] <- paste0("PC", 1:(ncol(cerw_dp15_maxmiss100_ld01_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp15_maxmiss100_ld01_loc <- rep(NA, length(cerw_dp15_maxmiss100_ld01_pca$ind))
cerw_dp15_maxmiss100_ld01_loc[grep("IN", cerw_dp15_maxmiss100_ld01_pca$ind)] <- "Indiana"
cerw_dp15_maxmiss100_ld01_loc[grep("TN", cerw_dp15_maxmiss100_ld01_pca$ind)] <- "Tennessee"
cerw_dp15_maxmiss100_ld01_loc[grep("OZ", cerw_dp15_maxmiss100_ld01_pca$ind)] <- "Ozarks"
cerw_dp15_maxmiss100_ld01_loc[grep("WI", cerw_dp15_maxmiss100_ld01_pca$ind)] <- "Wisconsin"
cerw_dp15_maxmiss100_ld01_loc[grep("PA", cerw_dp15_maxmiss100_ld01_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp15_maxmiss100_ld01_pca <- as_tibble(data.frame(cerw_dp15_maxmiss100_ld01_pca, cerw_dp15_maxmiss100_ld01_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp15_maxmiss100_ld01_pve <- data.frame(PC = 1:20, pve = cerw_dp15_maxmiss100_ld01_eigenval/sum(cerw_dp15_maxmiss100_ld01_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp15_maxmiss100_ld01_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp15_maxmiss100_ld01_pve$pve)

# plot pca
cerw_dp15_maxmiss100_ld01_pca_plot<-ggplot(cerw_dp15_maxmiss100_ld01_pca, aes(PC1, PC2, col = cerw_dp15_maxmiss100_ld01_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.1)) +
    xlab(paste0("PC1 (", signif(cerw_dp15_maxmiss100_ld01_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp15_maxmiss100_ld01_pve$pve[2], 3), "%)"))

cerw_dp15_maxmiss100_ld01_pca_plot



#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld02 ####
cerw_dp15_maxmiss100_ld02_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld02.eigenvec", col_names = FALSE)
cerw_dp15_maxmiss100_ld02_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld02.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp15_maxmiss100_ld02_pca <- cerw_dp15_maxmiss100_ld02_pca[,-1]
# set names
names(cerw_dp15_maxmiss100_ld02_pca)[1] <- "ind"
names(cerw_dp15_maxmiss100_ld02_pca)[2:ncol(cerw_dp15_maxmiss100_ld02_pca)] <- paste0("PC", 1:(ncol(cerw_dp15_maxmiss100_ld02_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp15_maxmiss100_ld02_loc <- rep(NA, length(cerw_dp15_maxmiss100_ld02_pca$ind))
cerw_dp15_maxmiss100_ld02_loc[grep("IN", cerw_dp15_maxmiss100_ld02_pca$ind)] <- "Indiana"
cerw_dp15_maxmiss100_ld02_loc[grep("TN", cerw_dp15_maxmiss100_ld02_pca$ind)] <- "Tennessee"
cerw_dp15_maxmiss100_ld02_loc[grep("OZ", cerw_dp15_maxmiss100_ld02_pca$ind)] <- "Ozarks"
cerw_dp15_maxmiss100_ld02_loc[grep("WI", cerw_dp15_maxmiss100_ld02_pca$ind)] <- "Wisconsin"
cerw_dp15_maxmiss100_ld02_loc[grep("PA", cerw_dp15_maxmiss100_ld02_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp15_maxmiss100_ld02_pca <- as_tibble(data.frame(cerw_dp15_maxmiss100_ld02_pca, cerw_dp15_maxmiss100_ld02_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp15_maxmiss100_ld02_pve <- data.frame(PC = 1:20, pve = cerw_dp15_maxmiss100_ld02_eigenval/sum(cerw_dp15_maxmiss100_ld02_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp15_maxmiss100_ld02_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp15_maxmiss100_ld02_pve$pve)

# plot pca
cerw_dp15_maxmiss100_ld02_pca_plot<-ggplot(cerw_dp15_maxmiss100_ld02_pca, aes(PC1, PC2, col = cerw_dp15_maxmiss100_ld02_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.2)) +
    xlab(paste0("PC1 (", signif(cerw_dp15_maxmiss100_ld02_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp15_maxmiss100_ld02_pve$pve[2], 3), "%)"))

cerw_dp15_maxmiss100_ld02_pca_plot


















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld03 ####
cerw_dp15_maxmiss100_ld03_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld03.eigenvec", col_names = FALSE)
cerw_dp15_maxmiss100_ld03_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld03.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp15_maxmiss100_ld03_pca <- cerw_dp15_maxmiss100_ld03_pca[,-1]
# set names
names(cerw_dp15_maxmiss100_ld03_pca)[1] <- "ind"
names(cerw_dp15_maxmiss100_ld03_pca)[2:ncol(cerw_dp15_maxmiss100_ld03_pca)] <- paste0("PC", 1:(ncol(cerw_dp15_maxmiss100_ld03_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp15_maxmiss100_ld03_loc <- rep(NA, length(cerw_dp15_maxmiss100_ld03_pca$ind))
cerw_dp15_maxmiss100_ld03_loc[grep("IN", cerw_dp15_maxmiss100_ld03_pca$ind)] <- "Indiana"
cerw_dp15_maxmiss100_ld03_loc[grep("TN", cerw_dp15_maxmiss100_ld03_pca$ind)] <- "Tennessee"
cerw_dp15_maxmiss100_ld03_loc[grep("OZ", cerw_dp15_maxmiss100_ld03_pca$ind)] <- "Ozarks"
cerw_dp15_maxmiss100_ld03_loc[grep("WI", cerw_dp15_maxmiss100_ld03_pca$ind)] <- "Wisconsin"
cerw_dp15_maxmiss100_ld03_loc[grep("PA", cerw_dp15_maxmiss100_ld03_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp15_maxmiss100_ld03_pca <- as_tibble(data.frame(cerw_dp15_maxmiss100_ld03_pca, cerw_dp15_maxmiss100_ld03_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp15_maxmiss100_ld03_pve <- data.frame(PC = 1:20, pve = cerw_dp15_maxmiss100_ld03_eigenval/sum(cerw_dp15_maxmiss100_ld03_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp15_maxmiss100_ld03_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp15_maxmiss100_ld03_pve$pve)

# plot pca
cerw_dp15_maxmiss100_ld03_pca_plot<-ggplot(cerw_dp15_maxmiss100_ld03_pca, aes(PC1, PC2, col = cerw_dp15_maxmiss100_ld03_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.3)) +
    xlab(paste0("PC1 (", signif(cerw_dp15_maxmiss100_ld03_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp15_maxmiss100_ld03_pve$pve[2], 3), "%)"))

cerw_dp15_maxmiss100_ld03_pca_plot















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld04 ####
cerw_dp15_maxmiss100_ld04_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld04.eigenvec", col_names = FALSE)
cerw_dp15_maxmiss100_ld04_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld04.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp15_maxmiss100_ld04_pca <- cerw_dp15_maxmiss100_ld04_pca[,-1]
# set names
names(cerw_dp15_maxmiss100_ld04_pca)[1] <- "ind"
names(cerw_dp15_maxmiss100_ld04_pca)[2:ncol(cerw_dp15_maxmiss100_ld04_pca)] <- paste0("PC", 1:(ncol(cerw_dp15_maxmiss100_ld04_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp15_maxmiss100_ld04_loc <- rep(NA, length(cerw_dp15_maxmiss100_ld04_pca$ind))
cerw_dp15_maxmiss100_ld04_loc[grep("IN", cerw_dp15_maxmiss100_ld04_pca$ind)] <- "Indiana"
cerw_dp15_maxmiss100_ld04_loc[grep("TN", cerw_dp15_maxmiss100_ld04_pca$ind)] <- "Tennessee"
cerw_dp15_maxmiss100_ld04_loc[grep("OZ", cerw_dp15_maxmiss100_ld04_pca$ind)] <- "Ozarks"
cerw_dp15_maxmiss100_ld04_loc[grep("WI", cerw_dp15_maxmiss100_ld04_pca$ind)] <- "Wisconsin"
cerw_dp15_maxmiss100_ld04_loc[grep("PA", cerw_dp15_maxmiss100_ld04_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp15_maxmiss100_ld04_pca <- as_tibble(data.frame(cerw_dp15_maxmiss100_ld04_pca, cerw_dp15_maxmiss100_ld04_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp15_maxmiss100_ld04_pve <- data.frame(PC = 1:20, pve = cerw_dp15_maxmiss100_ld04_eigenval/sum(cerw_dp15_maxmiss100_ld04_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp15_maxmiss100_ld04_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp15_maxmiss100_ld04_pve$pve)

# plot pca
cerw_dp15_maxmiss100_ld04_pca_plot<-ggplot(cerw_dp15_maxmiss100_ld04_pca, aes(PC1, PC2, col = cerw_dp15_maxmiss100_ld04_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.4)) +
    xlab(paste0("PC1 (", signif(cerw_dp15_maxmiss100_ld04_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp15_maxmiss100_ld04_pve$pve[2], 3), "%)"))

cerw_dp15_maxmiss100_ld04_pca_plot




#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld05 ####
cerw_dp15_maxmiss100_ld05_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld05.eigenvec", col_names = FALSE)
cerw_dp15_maxmiss100_ld05_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp15_maxmiss100_ld05_pca <- cerw_dp15_maxmiss100_ld05_pca[,-1]
# set names
names(cerw_dp15_maxmiss100_ld05_pca)[1] <- "ind"
names(cerw_dp15_maxmiss100_ld05_pca)[2:ncol(cerw_dp15_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(cerw_dp15_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp15_maxmiss100_ld05_loc <- rep(NA, length(cerw_dp15_maxmiss100_ld05_pca$ind))
cerw_dp15_maxmiss100_ld05_loc[grep("IN", cerw_dp15_maxmiss100_ld05_pca$ind)] <- "Indiana"
cerw_dp15_maxmiss100_ld05_loc[grep("TN", cerw_dp15_maxmiss100_ld05_pca$ind)] <- "Tennessee"
cerw_dp15_maxmiss100_ld05_loc[grep("OZ", cerw_dp15_maxmiss100_ld05_pca$ind)] <- "Ozarks"
cerw_dp15_maxmiss100_ld05_loc[grep("WI", cerw_dp15_maxmiss100_ld05_pca$ind)] <- "Wisconsin"
cerw_dp15_maxmiss100_ld05_loc[grep("PA", cerw_dp15_maxmiss100_ld05_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp15_maxmiss100_ld05_pca <- as_tibble(data.frame(cerw_dp15_maxmiss100_ld05_pca, cerw_dp15_maxmiss100_ld05_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp15_maxmiss100_ld05_pve <- data.frame(PC = 1:20, pve = cerw_dp15_maxmiss100_ld05_eigenval/sum(cerw_dp15_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp15_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp15_maxmiss100_ld05_pve$pve)

# plot pca
cerw_dp15_maxmiss100_ld05_pca_plot<-ggplot(cerw_dp15_maxmiss100_ld05_pca, aes(PC1, PC2, col = cerw_dp15_maxmiss100_ld05_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(cerw_dp15_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp15_maxmiss100_ld05_pve$pve[2], 3), "%)"))

cerw_dp15_maxmiss100_ld05_pca_plot











#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld08 ####
cerw_dp15_maxmiss100_ld08_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld08.eigenvec", col_names = FALSE)
cerw_dp15_maxmiss100_ld08_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld08.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp15_maxmiss100_ld08_pca <- cerw_dp15_maxmiss100_ld08_pca[,-1]
# set names
names(cerw_dp15_maxmiss100_ld08_pca)[1] <- "ind"
names(cerw_dp15_maxmiss100_ld08_pca)[2:ncol(cerw_dp15_maxmiss100_ld08_pca)] <- paste0("PC", 1:(ncol(cerw_dp15_maxmiss100_ld08_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp15_maxmiss100_ld08_loc <- rep(NA, length(cerw_dp15_maxmiss100_ld08_pca$ind))
cerw_dp15_maxmiss100_ld08_loc[grep("IN", cerw_dp15_maxmiss100_ld08_pca$ind)] <- "Indiana"
cerw_dp15_maxmiss100_ld08_loc[grep("TN", cerw_dp15_maxmiss100_ld08_pca$ind)] <- "Tennessee"
cerw_dp15_maxmiss100_ld08_loc[grep("OZ", cerw_dp15_maxmiss100_ld08_pca$ind)] <- "Ozarks"
cerw_dp15_maxmiss100_ld08_loc[grep("WI", cerw_dp15_maxmiss100_ld08_pca$ind)] <- "Wisconsin"
cerw_dp15_maxmiss100_ld08_loc[grep("PA", cerw_dp15_maxmiss100_ld08_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp15_maxmiss100_ld08_pca <- as_tibble(data.frame(cerw_dp15_maxmiss100_ld08_pca, cerw_dp15_maxmiss100_ld08_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp15_maxmiss100_ld08_pve <- data.frame(PC = 1:20, pve = cerw_dp15_maxmiss100_ld08_eigenval/sum(cerw_dp15_maxmiss100_ld08_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp15_maxmiss100_ld08_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp15_maxmiss100_ld08_pve$pve)

# plot pca
cerw_dp15_maxmiss100_ld08_pca_plot<-ggplot(cerw_dp15_maxmiss100_ld08_pca, aes(PC1, PC2, col = cerw_dp15_maxmiss100_ld08_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.8)) +
    xlab(paste0("PC1 (", signif(cerw_dp15_maxmiss100_ld08_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp15_maxmiss100_ld08_pve$pve[2], 3), "%)"))

cerw_dp15_maxmiss100_ld08_pca_plot



















#### CERW DEPTH SET AT 20X ####


#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld0001 ####
cerw_dp20_maxmiss100_ld0001_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld0001.eigenvec", col_names = FALSE)
cerw_dp20_maxmiss100_ld0001_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld0001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp20_maxmiss100_ld0001_pca <- cerw_dp20_maxmiss100_ld0001_pca[,-1]
# set names
names(cerw_dp20_maxmiss100_ld0001_pca)[1] <- "ind"

names(cerw_dp20_maxmiss100_ld0001_pca)[2:ncol(cerw_dp20_maxmiss100_ld0001_pca)] <- paste0("PC", 1:(ncol(cerw_dp20_maxmiss100_ld0001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp20_maxmiss100_ld0001_loc <- rep(NA, length(cerw_dp20_maxmiss100_ld0001_pca$ind))
cerw_dp20_maxmiss100_ld0001_loc[grep("IN", cerw_dp20_maxmiss100_ld0001_pca$ind)] <- "Indiana"
cerw_dp20_maxmiss100_ld0001_loc[grep("TN", cerw_dp20_maxmiss100_ld0001_pca$ind)] <- "Tennessee"
cerw_dp20_maxmiss100_ld0001_loc[grep("OZ", cerw_dp20_maxmiss100_ld0001_pca$ind)] <- "Ozarks"
cerw_dp20_maxmiss100_ld0001_loc[grep("WI", cerw_dp20_maxmiss100_ld0001_pca$ind)] <- "Wisconsin"
cerw_dp20_maxmiss100_ld0001_loc[grep("PA", cerw_dp20_maxmiss100_ld0001_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp20_maxmiss100_ld0001_pca <- as_tibble(data.frame(cerw_dp20_maxmiss100_ld0001_pca, cerw_dp20_maxmiss100_ld0001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp20_maxmiss100_ld0001_pve <- data.frame(PC = 1:20, pve = cerw_dp20_maxmiss100_ld0001_eigenval/sum(cerw_dp20_maxmiss100_ld0001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp20_maxmiss100_ld0001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp20_maxmiss100_ld0001_pve$pve)

# plot pca
cerw_dp20_maxmiss100_ld0001_pca_plot<-ggplot(cerw_dp20_maxmiss100_ld0001_pca, aes(PC1, PC2, col = cerw_dp20_maxmiss100_ld0001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.001)) +
    xlab(paste0("PC1 (", signif(cerw_dp20_maxmiss100_ld0001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp20_maxmiss100_ld0001_pve$pve[2], 3), "%)"))

cerw_dp20_maxmiss100_ld0001_pca_plot











#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld001 ####
cerw_dp20_maxmiss100_ld001_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld001.eigenvec", col_names = FALSE)
cerw_dp20_maxmiss100_ld001_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp20_maxmiss100_ld001_pca <- cerw_dp20_maxmiss100_ld001_pca[,-1]
# set names
names(cerw_dp20_maxmiss100_ld001_pca)[1] <- "ind"
names(cerw_dp20_maxmiss100_ld001_pca)[2:ncol(cerw_dp20_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(cerw_dp20_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp20_maxmiss100_ld001_loc <- rep(NA, length(cerw_dp20_maxmiss100_ld001_pca$ind))
cerw_dp20_maxmiss100_ld001_loc[grep("IN", cerw_dp20_maxmiss100_ld001_pca$ind)] <- "Indiana"
cerw_dp20_maxmiss100_ld001_loc[grep("TN", cerw_dp20_maxmiss100_ld001_pca$ind)] <- "Tennessee"
cerw_dp20_maxmiss100_ld001_loc[grep("OZ", cerw_dp20_maxmiss100_ld001_pca$ind)] <- "Ozarks"
cerw_dp20_maxmiss100_ld001_loc[grep("WI", cerw_dp20_maxmiss100_ld001_pca$ind)] <- "Wisconsin"
cerw_dp20_maxmiss100_ld001_loc[grep("PA", cerw_dp20_maxmiss100_ld001_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp20_maxmiss100_ld001_pca <- as_tibble(data.frame(cerw_dp20_maxmiss100_ld001_pca, cerw_dp20_maxmiss100_ld001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp20_maxmiss100_ld001_pve <- data.frame(PC = 1:20, pve = cerw_dp20_maxmiss100_ld001_eigenval/sum(cerw_dp20_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp20_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp20_maxmiss100_ld001_pve$pve)

# plot pca
cerw_dp20_maxmiss100_ld001_pca_plot<-ggplot(cerw_dp20_maxmiss100_ld001_pca, aes(PC1, PC2, col = cerw_dp20_maxmiss100_ld001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(cerw_dp20_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp20_maxmiss100_ld001_pve$pve[2], 3), "%)"))

cerw_dp20_maxmiss100_ld001_pca_plot















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld01 ####
cerw_dp20_maxmiss100_ld01_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld01.eigenvec", col_names = FALSE)
cerw_dp20_maxmiss100_ld01_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld01.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp20_maxmiss100_ld01_pca <- cerw_dp20_maxmiss100_ld01_pca[,-1]
# set names
names(cerw_dp20_maxmiss100_ld01_pca)[1] <- "ind"
names(cerw_dp20_maxmiss100_ld01_pca)[2:ncol(cerw_dp20_maxmiss100_ld01_pca)] <- paste0("PC", 1:(ncol(cerw_dp20_maxmiss100_ld01_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp20_maxmiss100_ld01_loc <- rep(NA, length(cerw_dp20_maxmiss100_ld01_pca$ind))
cerw_dp20_maxmiss100_ld01_loc[grep("IN", cerw_dp20_maxmiss100_ld01_pca$ind)] <- "Indiana"
cerw_dp20_maxmiss100_ld01_loc[grep("TN", cerw_dp20_maxmiss100_ld01_pca$ind)] <- "Tennessee"
cerw_dp20_maxmiss100_ld01_loc[grep("OZ", cerw_dp20_maxmiss100_ld01_pca$ind)] <- "Ozarks"
cerw_dp20_maxmiss100_ld01_loc[grep("WI", cerw_dp20_maxmiss100_ld01_pca$ind)] <- "Wisconsin"
cerw_dp20_maxmiss100_ld01_loc[grep("PA", cerw_dp20_maxmiss100_ld01_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp20_maxmiss100_ld01_pca <- as_tibble(data.frame(cerw_dp20_maxmiss100_ld01_pca, cerw_dp20_maxmiss100_ld01_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp20_maxmiss100_ld01_pve <- data.frame(PC = 1:20, pve = cerw_dp20_maxmiss100_ld01_eigenval/sum(cerw_dp20_maxmiss100_ld01_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp20_maxmiss100_ld01_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp20_maxmiss100_ld01_pve$pve)

# plot pca
cerw_dp20_maxmiss100_ld01_pca_plot<-ggplot(cerw_dp20_maxmiss100_ld01_pca, aes(PC1, PC2, col = cerw_dp20_maxmiss100_ld01_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.1)) +
    xlab(paste0("PC1 (", signif(cerw_dp20_maxmiss100_ld01_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp20_maxmiss100_ld01_pve$pve[2], 3), "%)"))

cerw_dp20_maxmiss100_ld01_pca_plot



#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld02 ####
cerw_dp20_maxmiss100_ld02_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld02.eigenvec", col_names = FALSE)
cerw_dp20_maxmiss100_ld02_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld02.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp20_maxmiss100_ld02_pca <- cerw_dp20_maxmiss100_ld02_pca[,-1]
# set names
names(cerw_dp20_maxmiss100_ld02_pca)[1] <- "ind"
names(cerw_dp20_maxmiss100_ld02_pca)[2:ncol(cerw_dp20_maxmiss100_ld02_pca)] <- paste0("PC", 1:(ncol(cerw_dp20_maxmiss100_ld02_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp20_maxmiss100_ld02_loc <- rep(NA, length(cerw_dp20_maxmiss100_ld02_pca$ind))
cerw_dp20_maxmiss100_ld02_loc[grep("IN", cerw_dp20_maxmiss100_ld02_pca$ind)] <- "Indiana"
cerw_dp20_maxmiss100_ld02_loc[grep("TN", cerw_dp20_maxmiss100_ld02_pca$ind)] <- "Tennessee"
cerw_dp20_maxmiss100_ld02_loc[grep("OZ", cerw_dp20_maxmiss100_ld02_pca$ind)] <- "Ozarks"
cerw_dp20_maxmiss100_ld02_loc[grep("WI", cerw_dp20_maxmiss100_ld02_pca$ind)] <- "Wisconsin"
cerw_dp20_maxmiss100_ld02_loc[grep("PA", cerw_dp20_maxmiss100_ld02_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp20_maxmiss100_ld02_pca <- as_tibble(data.frame(cerw_dp20_maxmiss100_ld02_pca, cerw_dp20_maxmiss100_ld02_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp20_maxmiss100_ld02_pve <- data.frame(PC = 1:20, pve = cerw_dp20_maxmiss100_ld02_eigenval/sum(cerw_dp20_maxmiss100_ld02_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp20_maxmiss100_ld02_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp20_maxmiss100_ld02_pve$pve)

# plot pca
cerw_dp20_maxmiss100_ld02_pca_plot<-ggplot(cerw_dp20_maxmiss100_ld02_pca, aes(PC1, PC2, col = cerw_dp20_maxmiss100_ld02_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.2)) +
    xlab(paste0("PC1 (", signif(cerw_dp20_maxmiss100_ld02_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp20_maxmiss100_ld02_pve$pve[2], 3), "%)"))

cerw_dp20_maxmiss100_ld02_pca_plot


















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld03 ####
cerw_dp20_maxmiss100_ld03_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld03.eigenvec", col_names = FALSE)
cerw_dp20_maxmiss100_ld03_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld03.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp20_maxmiss100_ld03_pca <- cerw_dp20_maxmiss100_ld03_pca[,-1]
# set names
names(cerw_dp20_maxmiss100_ld03_pca)[1] <- "ind"
names(cerw_dp20_maxmiss100_ld03_pca)[2:ncol(cerw_dp20_maxmiss100_ld03_pca)] <- paste0("PC", 1:(ncol(cerw_dp20_maxmiss100_ld03_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp20_maxmiss100_ld03_loc <- rep(NA, length(cerw_dp20_maxmiss100_ld03_pca$ind))
cerw_dp20_maxmiss100_ld03_loc[grep("IN", cerw_dp20_maxmiss100_ld03_pca$ind)] <- "Indiana"
cerw_dp20_maxmiss100_ld03_loc[grep("TN", cerw_dp20_maxmiss100_ld03_pca$ind)] <- "Tennessee"
cerw_dp20_maxmiss100_ld03_loc[grep("OZ", cerw_dp20_maxmiss100_ld03_pca$ind)] <- "Ozarks"
cerw_dp20_maxmiss100_ld03_loc[grep("WI", cerw_dp20_maxmiss100_ld03_pca$ind)] <- "Wisconsin"
cerw_dp20_maxmiss100_ld03_loc[grep("PA", cerw_dp20_maxmiss100_ld03_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp20_maxmiss100_ld03_pca <- as_tibble(data.frame(cerw_dp20_maxmiss100_ld03_pca, cerw_dp20_maxmiss100_ld03_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp20_maxmiss100_ld03_pve <- data.frame(PC = 1:20, pve = cerw_dp20_maxmiss100_ld03_eigenval/sum(cerw_dp20_maxmiss100_ld03_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp20_maxmiss100_ld03_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp20_maxmiss100_ld03_pve$pve)

# plot pca
cerw_dp20_maxmiss100_ld03_pca_plot<-ggplot(cerw_dp20_maxmiss100_ld03_pca, aes(PC1, PC2, col = cerw_dp20_maxmiss100_ld03_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.3)) +
    xlab(paste0("PC1 (", signif(cerw_dp20_maxmiss100_ld03_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp20_maxmiss100_ld03_pve$pve[2], 3), "%)"))

cerw_dp20_maxmiss100_ld03_pca_plot















#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld04 ####
cerw_dp20_maxmiss100_ld04_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld04.eigenvec", col_names = FALSE)
cerw_dp20_maxmiss100_ld04_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld04.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp20_maxmiss100_ld04_pca <- cerw_dp20_maxmiss100_ld04_pca[,-1]
# set names
names(cerw_dp20_maxmiss100_ld04_pca)[1] <- "ind"
names(cerw_dp20_maxmiss100_ld04_pca)[2:ncol(cerw_dp20_maxmiss100_ld04_pca)] <- paste0("PC", 1:(ncol(cerw_dp20_maxmiss100_ld04_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp20_maxmiss100_ld04_loc <- rep(NA, length(cerw_dp20_maxmiss100_ld04_pca$ind))
cerw_dp20_maxmiss100_ld04_loc[grep("IN", cerw_dp20_maxmiss100_ld04_pca$ind)] <- "Indiana"
cerw_dp20_maxmiss100_ld04_loc[grep("TN", cerw_dp20_maxmiss100_ld04_pca$ind)] <- "Tennessee"
cerw_dp20_maxmiss100_ld04_loc[grep("OZ", cerw_dp20_maxmiss100_ld04_pca$ind)] <- "Ozarks"
cerw_dp20_maxmiss100_ld04_loc[grep("WI", cerw_dp20_maxmiss100_ld04_pca$ind)] <- "Wisconsin"
cerw_dp20_maxmiss100_ld04_loc[grep("PA", cerw_dp20_maxmiss100_ld04_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp20_maxmiss100_ld04_pca <- as_tibble(data.frame(cerw_dp20_maxmiss100_ld04_pca, cerw_dp20_maxmiss100_ld04_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp20_maxmiss100_ld04_pve <- data.frame(PC = 1:20, pve = cerw_dp20_maxmiss100_ld04_eigenval/sum(cerw_dp20_maxmiss100_ld04_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp20_maxmiss100_ld04_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp20_maxmiss100_ld04_pve$pve)

# plot pca
cerw_dp20_maxmiss100_ld04_pca_plot<-ggplot(cerw_dp20_maxmiss100_ld04_pca, aes(PC1, PC2, col = cerw_dp20_maxmiss100_ld04_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.4)) +
    xlab(paste0("PC1 (", signif(cerw_dp20_maxmiss100_ld04_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp20_maxmiss100_ld04_pve$pve[2], 3), "%)"))

cerw_dp20_maxmiss100_ld04_pca_plot




#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld05 ####
cerw_dp20_maxmiss100_ld05_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld05.eigenvec", col_names = FALSE)
cerw_dp20_maxmiss100_ld05_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp20_maxmiss100_ld05_pca <- cerw_dp20_maxmiss100_ld05_pca[,-1]
# set names
names(cerw_dp20_maxmiss100_ld05_pca)[1] <- "ind"
names(cerw_dp20_maxmiss100_ld05_pca)[2:ncol(cerw_dp20_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(cerw_dp20_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp20_maxmiss100_ld05_loc <- rep(NA, length(cerw_dp20_maxmiss100_ld05_pca$ind))
cerw_dp20_maxmiss100_ld05_loc[grep("IN", cerw_dp20_maxmiss100_ld05_pca$ind)] <- "Indiana"
cerw_dp20_maxmiss100_ld05_loc[grep("TN", cerw_dp20_maxmiss100_ld05_pca$ind)] <- "Tennessee"
cerw_dp20_maxmiss100_ld05_loc[grep("OZ", cerw_dp20_maxmiss100_ld05_pca$ind)] <- "Ozarks"
cerw_dp20_maxmiss100_ld05_loc[grep("WI", cerw_dp20_maxmiss100_ld05_pca$ind)] <- "Wisconsin"
cerw_dp20_maxmiss100_ld05_loc[grep("PA", cerw_dp20_maxmiss100_ld05_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp20_maxmiss100_ld05_pca <- as_tibble(data.frame(cerw_dp20_maxmiss100_ld05_pca, cerw_dp20_maxmiss100_ld05_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp20_maxmiss100_ld05_pve <- data.frame(PC = 1:20, pve = cerw_dp20_maxmiss100_ld05_eigenval/sum(cerw_dp20_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp20_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp20_maxmiss100_ld05_pve$pve)

# plot pca
cerw_dp20_maxmiss100_ld05_pca_plot<-ggplot(cerw_dp20_maxmiss100_ld05_pca, aes(PC1, PC2, col = cerw_dp20_maxmiss100_ld05_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(cerw_dp20_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp20_maxmiss100_ld05_pve$pve[2], 3), "%)"))

cerw_dp20_maxmiss100_ld05_pca_plot








#### ++ CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld08 ####
cerw_dp20_maxmiss100_ld08_pca <- read_table("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld08.eigenvec", col_names = FALSE)
cerw_dp20_maxmiss100_ld08_eigenval <- scan("./inputs/CERW_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld08.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
cerw_dp20_maxmiss100_ld08_pca <- cerw_dp20_maxmiss100_ld08_pca[,-1]
# set names
names(cerw_dp20_maxmiss100_ld08_pca)[1] <- "ind"
names(cerw_dp20_maxmiss100_ld08_pca)[2:ncol(cerw_dp20_maxmiss100_ld08_pca)] <- paste0("PC", 1:(ncol(cerw_dp20_maxmiss100_ld08_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
cerw_dp20_maxmiss100_ld08_loc <- rep(NA, length(cerw_dp20_maxmiss100_ld08_pca$ind))
cerw_dp20_maxmiss100_ld08_loc[grep("IN", cerw_dp20_maxmiss100_ld08_pca$ind)] <- "Indiana"
cerw_dp20_maxmiss100_ld08_loc[grep("TN", cerw_dp20_maxmiss100_ld08_pca$ind)] <- "Tennessee"
cerw_dp20_maxmiss100_ld08_loc[grep("OZ", cerw_dp20_maxmiss100_ld08_pca$ind)] <- "Ozarks"
cerw_dp20_maxmiss100_ld08_loc[grep("WI", cerw_dp20_maxmiss100_ld08_pca$ind)] <- "Wisconsin"
cerw_dp20_maxmiss100_ld08_loc[grep("PA", cerw_dp20_maxmiss100_ld08_pca$ind)] <- "Pennsylvania"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
cerw_dp20_maxmiss100_ld08_pca <- as_tibble(data.frame(cerw_dp20_maxmiss100_ld08_pca, cerw_dp20_maxmiss100_ld08_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
cerw_dp20_maxmiss100_ld08_pve <- data.frame(PC = 1:20, pve = cerw_dp20_maxmiss100_ld08_eigenval/sum(cerw_dp20_maxmiss100_ld08_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(cerw_dp20_maxmiss100_ld08_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(cerw_dp20_maxmiss100_ld08_pve$pve)

# plot pca
cerw_dp20_maxmiss100_ld08_pca_plot<-ggplot(cerw_dp20_maxmiss100_ld08_pca, aes(PC1, PC2, col = cerw_dp20_maxmiss100_ld08_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. ischyros)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.8)) +
    xlab(paste0("PC1 (", signif(cerw_dp20_maxmiss100_ld08_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(cerw_dp20_maxmiss100_ld08_pve$pve[2], 3), "%)"))

cerw_dp20_maxmiss100_ld08_pca_plot














#### PROW DEPTH SET AT 5X ####



#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld0001 ####
prow_dp05_maxmiss100_ld0001_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld0001.eigenvec", col_names = FALSE)
prow_dp05_maxmiss100_ld0001_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld0001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp05_maxmiss100_ld0001_pca <- prow_dp05_maxmiss100_ld0001_pca[,-1]
# set names
names(prow_dp05_maxmiss100_ld0001_pca)[1] <- "ind"
names(prow_dp05_maxmiss100_ld0001_pca)[2:ncol(prow_dp05_maxmiss100_ld0001_pca)] <- paste0("PC", 1:(ncol(prow_dp05_maxmiss100_ld0001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp05_maxmiss100_ld0001_loc <- rep(NA, length(prow_dp05_maxmiss100_ld0001_pca$ind))
prow_dp05_maxmiss100_ld0001_loc[grep("VA", prow_dp05_maxmiss100_ld0001_pca$ind)] <- "Virginia"
prow_dp05_maxmiss100_ld0001_loc[grep("OH", prow_dp05_maxmiss100_ld0001_pca$ind)] <- "Ohio"
prow_dp05_maxmiss100_ld0001_loc[grep("AR", prow_dp05_maxmiss100_ld0001_pca$ind)] <- "Arkansas"
prow_dp05_maxmiss100_ld0001_loc[grep("WI", prow_dp05_maxmiss100_ld0001_pca$ind)] <- "Wisconsin"
prow_dp05_maxmiss100_ld0001_loc[grep("SC", prow_dp05_maxmiss100_ld0001_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp05_maxmiss100_ld0001_pca <- as_tibble(data.frame(prow_dp05_maxmiss100_ld0001_pca, prow_dp05_maxmiss100_ld0001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp05_maxmiss100_ld0001_pve <- data.frame(PC = 1:20, pve = prow_dp05_maxmiss100_ld0001_eigenval/sum(prow_dp05_maxmiss100_ld0001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp05_maxmiss100_ld0001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp05_maxmiss100_ld0001_pve$pve)

# plot pca
prow_dp05_maxmiss100_ld0001_pca_plot<-ggplot(prow_dp05_maxmiss100_ld0001_pca, aes(PC1, PC2, col = prow_dp05_maxmiss100_ld0001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.001)) +
    xlab(paste0("PC1 (", signif(prow_dp05_maxmiss100_ld0001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp05_maxmiss100_ld0001_pve$pve[2], 3), "%)"))

prow_dp05_maxmiss100_ld0001_pca_plot















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld001 ####
prow_dp05_maxmiss100_ld001_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_dp05_maxmiss100_ld001_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp05_maxmiss100_ld001_pca <- prow_dp05_maxmiss100_ld001_pca[,-1]
# set names
names(prow_dp05_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_dp05_maxmiss100_ld001_pca)[2:ncol(prow_dp05_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_dp05_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp05_maxmiss100_ld001_loc <- rep(NA, length(prow_dp05_maxmiss100_ld001_pca$ind))
prow_dp05_maxmiss100_ld001_loc[grep("VA", prow_dp05_maxmiss100_ld001_pca$ind)] <- "Virginia"
prow_dp05_maxmiss100_ld001_loc[grep("OH", prow_dp05_maxmiss100_ld001_pca$ind)] <- "Ohio"
prow_dp05_maxmiss100_ld001_loc[grep("AR", prow_dp05_maxmiss100_ld001_pca$ind)] <- "Arkansas"
prow_dp05_maxmiss100_ld001_loc[grep("WI", prow_dp05_maxmiss100_ld001_pca$ind)] <- "Wisconsin"
prow_dp05_maxmiss100_ld001_loc[grep("SC", prow_dp05_maxmiss100_ld001_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp05_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_dp05_maxmiss100_ld001_pca, prow_dp05_maxmiss100_ld001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp05_maxmiss100_ld001_pve <- data.frame(PC = 1:20, pve = prow_dp05_maxmiss100_ld001_eigenval/sum(prow_dp05_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp05_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp05_maxmiss100_ld001_pve$pve)

# plot pca
prow_dp05_maxmiss100_ld001_pca_plot<-ggplot(prow_dp05_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_dp05_maxmiss100_ld001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_dp05_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp05_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_dp05_maxmiss100_ld001_pca_plot




















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld01 ####
prow_dp05_maxmiss100_ld01_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld01.eigenvec", col_names = FALSE)
prow_dp05_maxmiss100_ld01_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld01.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp05_maxmiss100_ld01_pca <- prow_dp05_maxmiss100_ld01_pca[,-1]
# set names
names(prow_dp05_maxmiss100_ld01_pca)[1] <- "ind"
names(prow_dp05_maxmiss100_ld01_pca)[2:ncol(prow_dp05_maxmiss100_ld01_pca)] <- paste0("PC", 1:(ncol(prow_dp05_maxmiss100_ld01_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp05_maxmiss100_ld01_loc <- rep(NA, length(prow_dp05_maxmiss100_ld01_pca$ind))
prow_dp05_maxmiss100_ld01_loc[grep("VA", prow_dp05_maxmiss100_ld01_pca$ind)] <- "Virginia"
prow_dp05_maxmiss100_ld01_loc[grep("OH", prow_dp05_maxmiss100_ld01_pca$ind)] <- "Ohio"
prow_dp05_maxmiss100_ld01_loc[grep("AR", prow_dp05_maxmiss100_ld01_pca$ind)] <- "Arkansas"
prow_dp05_maxmiss100_ld01_loc[grep("WI", prow_dp05_maxmiss100_ld01_pca$ind)] <- "Wisconsin"
prow_dp05_maxmiss100_ld01_loc[grep("SC", prow_dp05_maxmiss100_ld01_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp05_maxmiss100_ld01_pca <- as_tibble(data.frame(prow_dp05_maxmiss100_ld01_pca, prow_dp05_maxmiss100_ld01_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp05_maxmiss100_ld01_pve <- data.frame(PC = 1:20, pve = prow_dp05_maxmiss100_ld01_eigenval/sum(prow_dp05_maxmiss100_ld01_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp05_maxmiss100_ld01_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp05_maxmiss100_ld01_pve$pve)

# plot pca
prow_dp05_maxmiss100_ld01_pca_plot<-ggplot(prow_dp05_maxmiss100_ld01_pca, aes(PC1, PC2, col = prow_dp05_maxmiss100_ld01_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.1)) +
    xlab(paste0("PC1 (", signif(prow_dp05_maxmiss100_ld01_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp05_maxmiss100_ld01_pve$pve[2], 3), "%)"))

prow_dp05_maxmiss100_ld01_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld02 ####
prow_dp05_maxmiss100_ld02_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld02.eigenvec", col_names = FALSE)
prow_dp05_maxmiss100_ld02_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld02.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp05_maxmiss100_ld02_pca <- prow_dp05_maxmiss100_ld02_pca[,-1]
# set names
names(prow_dp05_maxmiss100_ld02_pca)[1] <- "ind"
names(prow_dp05_maxmiss100_ld02_pca)[2:ncol(prow_dp05_maxmiss100_ld02_pca)] <- paste0("PC", 1:(ncol(prow_dp05_maxmiss100_ld02_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp05_maxmiss100_ld02_loc <- rep(NA, length(prow_dp05_maxmiss100_ld02_pca$ind))
prow_dp05_maxmiss100_ld02_loc[grep("VA", prow_dp05_maxmiss100_ld02_pca$ind)] <- "Virginia"
prow_dp05_maxmiss100_ld02_loc[grep("OH", prow_dp05_maxmiss100_ld02_pca$ind)] <- "Ohio"
prow_dp05_maxmiss100_ld02_loc[grep("AR", prow_dp05_maxmiss100_ld02_pca$ind)] <- "Arkansas"
prow_dp05_maxmiss100_ld02_loc[grep("WI", prow_dp05_maxmiss100_ld02_pca$ind)] <- "Wisconsin"
prow_dp05_maxmiss100_ld02_loc[grep("SC", prow_dp05_maxmiss100_ld02_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp05_maxmiss100_ld02_pca <- as_tibble(data.frame(prow_dp05_maxmiss100_ld02_pca, prow_dp05_maxmiss100_ld02_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp05_maxmiss100_ld02_pve <- data.frame(PC = 1:20, pve = prow_dp05_maxmiss100_ld02_eigenval/sum(prow_dp05_maxmiss100_ld02_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp05_maxmiss100_ld02_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp05_maxmiss100_ld02_pve$pve)

# plot pca
prow_dp05_maxmiss100_ld02_pca_plot<-ggplot(prow_dp05_maxmiss100_ld02_pca, aes(PC1, PC2, col = prow_dp05_maxmiss100_ld02_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.2)) +
    xlab(paste0("PC1 (", signif(prow_dp05_maxmiss100_ld02_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp05_maxmiss100_ld02_pve$pve[2], 3), "%)"))

prow_dp05_maxmiss100_ld02_pca_plot















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld03 ####
prow_dp05_maxmiss100_ld03_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld03.eigenvec", col_names = FALSE)
prow_dp05_maxmiss100_ld03_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld03.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp05_maxmiss100_ld03_pca <- prow_dp05_maxmiss100_ld03_pca[,-1]
# set names
names(prow_dp05_maxmiss100_ld03_pca)[1] <- "ind"
names(prow_dp05_maxmiss100_ld03_pca)[2:ncol(prow_dp05_maxmiss100_ld03_pca)] <- paste0("PC", 1:(ncol(prow_dp05_maxmiss100_ld03_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp05_maxmiss100_ld03_loc <- rep(NA, length(prow_dp05_maxmiss100_ld03_pca$ind))
prow_dp05_maxmiss100_ld03_loc[grep("VA", prow_dp05_maxmiss100_ld03_pca$ind)] <- "Virginia"
prow_dp05_maxmiss100_ld03_loc[grep("OH", prow_dp05_maxmiss100_ld03_pca$ind)] <- "Ohio"
prow_dp05_maxmiss100_ld03_loc[grep("AR", prow_dp05_maxmiss100_ld03_pca$ind)] <- "Arkansas"
prow_dp05_maxmiss100_ld03_loc[grep("WI", prow_dp05_maxmiss100_ld03_pca$ind)] <- "Wisconsin"
prow_dp05_maxmiss100_ld03_loc[grep("SC", prow_dp05_maxmiss100_ld03_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp05_maxmiss100_ld03_pca <- as_tibble(data.frame(prow_dp05_maxmiss100_ld03_pca, prow_dp05_maxmiss100_ld03_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp05_maxmiss100_ld03_pve <- data.frame(PC = 1:20, pve = prow_dp05_maxmiss100_ld03_eigenval/sum(prow_dp05_maxmiss100_ld03_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp05_maxmiss100_ld03_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp05_maxmiss100_ld03_pve$pve)

# plot pca
prow_dp05_maxmiss100_ld03_pca_plot<-ggplot(prow_dp05_maxmiss100_ld03_pca, aes(PC1, PC2, col = prow_dp05_maxmiss100_ld03_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.3)) +
    xlab(paste0("PC1 (", signif(prow_dp05_maxmiss100_ld03_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp05_maxmiss100_ld03_pve$pve[2], 3), "%)"))

prow_dp05_maxmiss100_ld03_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld04 ####
prow_dp05_maxmiss100_ld04_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld04.eigenvec", col_names = FALSE)
prow_dp05_maxmiss100_ld04_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld04.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp05_maxmiss100_ld04_pca <- prow_dp05_maxmiss100_ld04_pca[,-1]
# set names
names(prow_dp05_maxmiss100_ld04_pca)[1] <- "ind"
names(prow_dp05_maxmiss100_ld04_pca)[2:ncol(prow_dp05_maxmiss100_ld04_pca)] <- paste0("PC", 1:(ncol(prow_dp05_maxmiss100_ld04_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp05_maxmiss100_ld04_loc <- rep(NA, length(prow_dp05_maxmiss100_ld04_pca$ind))
prow_dp05_maxmiss100_ld04_loc[grep("VA", prow_dp05_maxmiss100_ld04_pca$ind)] <- "Virginia"
prow_dp05_maxmiss100_ld04_loc[grep("OH", prow_dp05_maxmiss100_ld04_pca$ind)] <- "Ohio"
prow_dp05_maxmiss100_ld04_loc[grep("AR", prow_dp05_maxmiss100_ld04_pca$ind)] <- "Arkansas"
prow_dp05_maxmiss100_ld04_loc[grep("WI", prow_dp05_maxmiss100_ld04_pca$ind)] <- "Wisconsin"
prow_dp05_maxmiss100_ld04_loc[grep("SC", prow_dp05_maxmiss100_ld04_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp05_maxmiss100_ld04_pca <- as_tibble(data.frame(prow_dp05_maxmiss100_ld04_pca, prow_dp05_maxmiss100_ld04_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp05_maxmiss100_ld04_pve <- data.frame(PC = 1:20, pve = prow_dp05_maxmiss100_ld04_eigenval/sum(prow_dp05_maxmiss100_ld04_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp05_maxmiss100_ld04_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp05_maxmiss100_ld04_pve$pve)

# plot pca
prow_dp05_maxmiss100_ld04_pca_plot<-ggplot(prow_dp05_maxmiss100_ld04_pca, aes(PC1, PC2, col = prow_dp05_maxmiss100_ld04_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.4)) +
    xlab(paste0("PC1 (", signif(prow_dp05_maxmiss100_ld04_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp05_maxmiss100_ld04_pve$pve[2], 3), "%)"))

prow_dp05_maxmiss100_ld04_pca_plot





















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld05 ####
prow_dp05_maxmiss100_ld05_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_dp05_maxmiss100_ld05_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp05_maxmiss100_ld05_pca <- prow_dp05_maxmiss100_ld05_pca[,-1]
# set names
names(prow_dp05_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_dp05_maxmiss100_ld05_pca)[2:ncol(prow_dp05_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_dp05_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp05_maxmiss100_ld05_loc <- rep(NA, length(prow_dp05_maxmiss100_ld05_pca$ind))
prow_dp05_maxmiss100_ld05_loc[grep("VA", prow_dp05_maxmiss100_ld05_pca$ind)] <- "Virginia"
prow_dp05_maxmiss100_ld05_loc[grep("OH", prow_dp05_maxmiss100_ld05_pca$ind)] <- "Ohio"
prow_dp05_maxmiss100_ld05_loc[grep("AR", prow_dp05_maxmiss100_ld05_pca$ind)] <- "Arkansas"
prow_dp05_maxmiss100_ld05_loc[grep("WI", prow_dp05_maxmiss100_ld05_pca$ind)] <- "Wisconsin"
prow_dp05_maxmiss100_ld05_loc[grep("SC", prow_dp05_maxmiss100_ld05_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp05_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_dp05_maxmiss100_ld05_pca, prow_dp05_maxmiss100_ld05_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp05_maxmiss100_ld05_pve <- data.frame(PC = 1:20, pve = prow_dp05_maxmiss100_ld05_eigenval/sum(prow_dp05_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp05_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp05_maxmiss100_ld05_pve$pve)

# plot pca
prow_dp05_maxmiss100_ld05_pca_plot<-ggplot(prow_dp05_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_dp05_maxmiss100_ld05_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_dp05_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp05_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_dp05_maxmiss100_ld05_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld08 ####
prow_dp05_maxmiss100_ld08_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld08.eigenvec", col_names = FALSE)
prow_dp05_maxmiss100_ld08_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp05_maxmiss100_ld08.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp05_maxmiss100_ld08_pca <- prow_dp05_maxmiss100_ld08_pca[,-1]
# set names
names(prow_dp05_maxmiss100_ld08_pca)[1] <- "ind"
names(prow_dp05_maxmiss100_ld08_pca)[2:ncol(prow_dp05_maxmiss100_ld08_pca)] <- paste0("PC", 1:(ncol(prow_dp05_maxmiss100_ld08_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp05_maxmiss100_ld08_loc <- rep(NA, length(prow_dp05_maxmiss100_ld08_pca$ind))
prow_dp05_maxmiss100_ld08_loc[grep("VA", prow_dp05_maxmiss100_ld08_pca$ind)] <- "Virginia"
prow_dp05_maxmiss100_ld08_loc[grep("OH", prow_dp05_maxmiss100_ld08_pca$ind)] <- "Ohio"
prow_dp05_maxmiss100_ld08_loc[grep("AR", prow_dp05_maxmiss100_ld08_pca$ind)] <- "Arkansas"
prow_dp05_maxmiss100_ld08_loc[grep("WI", prow_dp05_maxmiss100_ld08_pca$ind)] <- "Wisconsin"
prow_dp05_maxmiss100_ld08_loc[grep("SC", prow_dp05_maxmiss100_ld08_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp05_maxmiss100_ld08_pca <- as_tibble(data.frame(prow_dp05_maxmiss100_ld08_pca, prow_dp05_maxmiss100_ld08_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp05_maxmiss100_ld08_pve <- data.frame(PC = 1:20, pve = prow_dp05_maxmiss100_ld08_eigenval/sum(prow_dp05_maxmiss100_ld08_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp05_maxmiss100_ld08_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp05_maxmiss100_ld08_pve$pve)

# plot pca
prow_dp05_maxmiss100_ld08_pca_plot<-ggplot(prow_dp05_maxmiss100_ld08_pca, aes(PC1, PC2, col = prow_dp05_maxmiss100_ld08_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.8)) +
    xlab(paste0("PC1 (", signif(prow_dp05_maxmiss100_ld08_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp05_maxmiss100_ld08_pve$pve[2], 3), "%)"))

prow_dp05_maxmiss100_ld08_pca_plot























#### PROW DEPTH SET AT 10X ####


#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld0001 ####
prow_dp10_maxmiss100_ld0001_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld0001.eigenvec", col_names = FALSE)
prow_dp10_maxmiss100_ld0001_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld0001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp10_maxmiss100_ld0001_pca <- prow_dp10_maxmiss100_ld0001_pca[,-1]
# set names
names(prow_dp10_maxmiss100_ld0001_pca)[1] <- "ind"
names(prow_dp10_maxmiss100_ld0001_pca)[2:ncol(prow_dp10_maxmiss100_ld0001_pca)] <- paste0("PC", 1:(ncol(prow_dp10_maxmiss100_ld0001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp10_maxmiss100_ld0001_loc <- rep(NA, length(prow_dp10_maxmiss100_ld0001_pca$ind))
prow_dp10_maxmiss100_ld0001_loc[grep("VA", prow_dp10_maxmiss100_ld0001_pca$ind)] <- "Virginia"
prow_dp10_maxmiss100_ld0001_loc[grep("OH", prow_dp10_maxmiss100_ld0001_pca$ind)] <- "Ohio"
prow_dp10_maxmiss100_ld0001_loc[grep("AR", prow_dp10_maxmiss100_ld0001_pca$ind)] <- "Arkansas"
prow_dp10_maxmiss100_ld0001_loc[grep("WI", prow_dp10_maxmiss100_ld0001_pca$ind)] <- "Wisconsin"
prow_dp10_maxmiss100_ld0001_loc[grep("SC", prow_dp10_maxmiss100_ld0001_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp10_maxmiss100_ld0001_pca <- as_tibble(data.frame(prow_dp10_maxmiss100_ld0001_pca, prow_dp10_maxmiss100_ld0001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp10_maxmiss100_ld0001_pve <- data.frame(PC = 1:20, pve = prow_dp10_maxmiss100_ld0001_eigenval/sum(prow_dp10_maxmiss100_ld0001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp10_maxmiss100_ld0001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp10_maxmiss100_ld0001_pve$pve)

# plot pca
prow_dp10_maxmiss100_ld0001_pca_plot<-ggplot(prow_dp10_maxmiss100_ld0001_pca, aes(PC1, PC2, col = prow_dp10_maxmiss100_ld0001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.001)) +
    xlab(paste0("PC1 (", signif(prow_dp10_maxmiss100_ld0001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp10_maxmiss100_ld0001_pve$pve[2], 3), "%)"))

prow_dp10_maxmiss100_ld0001_pca_plot















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001 ####
prow_dp10_maxmiss100_ld001_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_dp10_maxmiss100_ld001_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp10_maxmiss100_ld001_pca <- prow_dp10_maxmiss100_ld001_pca[,-1]
# set names
names(prow_dp10_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_dp10_maxmiss100_ld001_pca)[2:ncol(prow_dp10_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_dp10_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp10_maxmiss100_ld001_loc <- rep(NA, length(prow_dp10_maxmiss100_ld001_pca$ind))
prow_dp10_maxmiss100_ld001_loc[grep("VA", prow_dp10_maxmiss100_ld001_pca$ind)] <- "Virginia"
prow_dp10_maxmiss100_ld001_loc[grep("OH", prow_dp10_maxmiss100_ld001_pca$ind)] <- "Ohio"
prow_dp10_maxmiss100_ld001_loc[grep("AR", prow_dp10_maxmiss100_ld001_pca$ind)] <- "Arkansas"
prow_dp10_maxmiss100_ld001_loc[grep("WI", prow_dp10_maxmiss100_ld001_pca$ind)] <- "Wisconsin"
prow_dp10_maxmiss100_ld001_loc[grep("SC", prow_dp10_maxmiss100_ld001_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp10_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_dp10_maxmiss100_ld001_pca, prow_dp10_maxmiss100_ld001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp10_maxmiss100_ld001_pve <- data.frame(PC = 1:20, pve = prow_dp10_maxmiss100_ld001_eigenval/sum(prow_dp10_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp10_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp10_maxmiss100_ld001_pve$pve)

# plot pca
prow_dp10_maxmiss100_ld001_pca_plot<-ggplot(prow_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_dp10_maxmiss100_ld001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_dp10_maxmiss100_ld001_pca_plot




















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld01 ####
prow_dp10_maxmiss100_ld01_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld01.eigenvec", col_names = FALSE)
prow_dp10_maxmiss100_ld01_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld01.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp10_maxmiss100_ld01_pca <- prow_dp10_maxmiss100_ld01_pca[,-1]
# set names
names(prow_dp10_maxmiss100_ld01_pca)[1] <- "ind"
names(prow_dp10_maxmiss100_ld01_pca)[2:ncol(prow_dp10_maxmiss100_ld01_pca)] <- paste0("PC", 1:(ncol(prow_dp10_maxmiss100_ld01_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp10_maxmiss100_ld01_loc <- rep(NA, length(prow_dp10_maxmiss100_ld01_pca$ind))
prow_dp10_maxmiss100_ld01_loc[grep("VA", prow_dp10_maxmiss100_ld01_pca$ind)] <- "Virginia"
prow_dp10_maxmiss100_ld01_loc[grep("OH", prow_dp10_maxmiss100_ld01_pca$ind)] <- "Ohio"
prow_dp10_maxmiss100_ld01_loc[grep("AR", prow_dp10_maxmiss100_ld01_pca$ind)] <- "Arkansas"
prow_dp10_maxmiss100_ld01_loc[grep("WI", prow_dp10_maxmiss100_ld01_pca$ind)] <- "Wisconsin"
prow_dp10_maxmiss100_ld01_loc[grep("SC", prow_dp10_maxmiss100_ld01_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp10_maxmiss100_ld01_pca <- as_tibble(data.frame(prow_dp10_maxmiss100_ld01_pca, prow_dp10_maxmiss100_ld01_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp10_maxmiss100_ld01_pve <- data.frame(PC = 1:20, pve = prow_dp10_maxmiss100_ld01_eigenval/sum(prow_dp10_maxmiss100_ld01_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp10_maxmiss100_ld01_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp10_maxmiss100_ld01_pve$pve)

# plot pca
prow_dp10_maxmiss100_ld01_pca_plot<-ggplot(prow_dp10_maxmiss100_ld01_pca, aes(PC1, PC2, col = prow_dp10_maxmiss100_ld01_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.1)) +
    xlab(paste0("PC1 (", signif(prow_dp10_maxmiss100_ld01_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp10_maxmiss100_ld01_pve$pve[2], 3), "%)"))

prow_dp10_maxmiss100_ld01_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld02 ####
prow_dp10_maxmiss100_ld02_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld02.eigenvec", col_names = FALSE)
prow_dp10_maxmiss100_ld02_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld02.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp10_maxmiss100_ld02_pca <- prow_dp10_maxmiss100_ld02_pca[,-1]
# set names
names(prow_dp10_maxmiss100_ld02_pca)[1] <- "ind"
names(prow_dp10_maxmiss100_ld02_pca)[2:ncol(prow_dp10_maxmiss100_ld02_pca)] <- paste0("PC", 1:(ncol(prow_dp10_maxmiss100_ld02_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp10_maxmiss100_ld02_loc <- rep(NA, length(prow_dp10_maxmiss100_ld02_pca$ind))
prow_dp10_maxmiss100_ld02_loc[grep("VA", prow_dp10_maxmiss100_ld02_pca$ind)] <- "Virginia"
prow_dp10_maxmiss100_ld02_loc[grep("OH", prow_dp10_maxmiss100_ld02_pca$ind)] <- "Ohio"
prow_dp10_maxmiss100_ld02_loc[grep("AR", prow_dp10_maxmiss100_ld02_pca$ind)] <- "Arkansas"
prow_dp10_maxmiss100_ld02_loc[grep("WI", prow_dp10_maxmiss100_ld02_pca$ind)] <- "Wisconsin"
prow_dp10_maxmiss100_ld02_loc[grep("SC", prow_dp10_maxmiss100_ld02_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp10_maxmiss100_ld02_pca <- as_tibble(data.frame(prow_dp10_maxmiss100_ld02_pca, prow_dp10_maxmiss100_ld02_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp10_maxmiss100_ld02_pve <- data.frame(PC = 1:20, pve = prow_dp10_maxmiss100_ld02_eigenval/sum(prow_dp10_maxmiss100_ld02_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp10_maxmiss100_ld02_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp10_maxmiss100_ld02_pve$pve)

# plot pca
prow_dp10_maxmiss100_ld02_pca_plot<-ggplot(prow_dp10_maxmiss100_ld02_pca, aes(PC1, PC2, col = prow_dp10_maxmiss100_ld02_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.2)) +
    xlab(paste0("PC1 (", signif(prow_dp10_maxmiss100_ld02_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp10_maxmiss100_ld02_pve$pve[2], 3), "%)"))

prow_dp10_maxmiss100_ld02_pca_plot















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld03 ####
prow_dp10_maxmiss100_ld03_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld03.eigenvec", col_names = FALSE)
prow_dp10_maxmiss100_ld03_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld03.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp10_maxmiss100_ld03_pca <- prow_dp10_maxmiss100_ld03_pca[,-1]
# set names
names(prow_dp10_maxmiss100_ld03_pca)[1] <- "ind"
names(prow_dp10_maxmiss100_ld03_pca)[2:ncol(prow_dp10_maxmiss100_ld03_pca)] <- paste0("PC", 1:(ncol(prow_dp10_maxmiss100_ld03_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp10_maxmiss100_ld03_loc <- rep(NA, length(prow_dp10_maxmiss100_ld03_pca$ind))
prow_dp10_maxmiss100_ld03_loc[grep("VA", prow_dp10_maxmiss100_ld03_pca$ind)] <- "Virginia"
prow_dp10_maxmiss100_ld03_loc[grep("OH", prow_dp10_maxmiss100_ld03_pca$ind)] <- "Ohio"
prow_dp10_maxmiss100_ld03_loc[grep("AR", prow_dp10_maxmiss100_ld03_pca$ind)] <- "Arkansas"
prow_dp10_maxmiss100_ld03_loc[grep("WI", prow_dp10_maxmiss100_ld03_pca$ind)] <- "Wisconsin"
prow_dp10_maxmiss100_ld03_loc[grep("SC", prow_dp10_maxmiss100_ld03_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp10_maxmiss100_ld03_pca <- as_tibble(data.frame(prow_dp10_maxmiss100_ld03_pca, prow_dp10_maxmiss100_ld03_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp10_maxmiss100_ld03_pve <- data.frame(PC = 1:20, pve = prow_dp10_maxmiss100_ld03_eigenval/sum(prow_dp10_maxmiss100_ld03_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp10_maxmiss100_ld03_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp10_maxmiss100_ld03_pve$pve)

# plot pca
prow_dp10_maxmiss100_ld03_pca_plot<-ggplot(prow_dp10_maxmiss100_ld03_pca, aes(PC1, PC2, col = prow_dp10_maxmiss100_ld03_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.3)) +
    xlab(paste0("PC1 (", signif(prow_dp10_maxmiss100_ld03_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp10_maxmiss100_ld03_pve$pve[2], 3), "%)"))

prow_dp10_maxmiss100_ld03_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld04 ####
prow_dp10_maxmiss100_ld04_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld04.eigenvec", col_names = FALSE)
prow_dp10_maxmiss100_ld04_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld04.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp10_maxmiss100_ld04_pca <- prow_dp10_maxmiss100_ld04_pca[,-1]
# set names
names(prow_dp10_maxmiss100_ld04_pca)[1] <- "ind"
names(prow_dp10_maxmiss100_ld04_pca)[2:ncol(prow_dp10_maxmiss100_ld04_pca)] <- paste0("PC", 1:(ncol(prow_dp10_maxmiss100_ld04_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp10_maxmiss100_ld04_loc <- rep(NA, length(prow_dp10_maxmiss100_ld04_pca$ind))
prow_dp10_maxmiss100_ld04_loc[grep("VA", prow_dp10_maxmiss100_ld04_pca$ind)] <- "Virginia"
prow_dp10_maxmiss100_ld04_loc[grep("OH", prow_dp10_maxmiss100_ld04_pca$ind)] <- "Ohio"
prow_dp10_maxmiss100_ld04_loc[grep("AR", prow_dp10_maxmiss100_ld04_pca$ind)] <- "Arkansas"
prow_dp10_maxmiss100_ld04_loc[grep("WI", prow_dp10_maxmiss100_ld04_pca$ind)] <- "Wisconsin"
prow_dp10_maxmiss100_ld04_loc[grep("SC", prow_dp10_maxmiss100_ld04_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp10_maxmiss100_ld04_pca <- as_tibble(data.frame(prow_dp10_maxmiss100_ld04_pca, prow_dp10_maxmiss100_ld04_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp10_maxmiss100_ld04_pve <- data.frame(PC = 1:20, pve = prow_dp10_maxmiss100_ld04_eigenval/sum(prow_dp10_maxmiss100_ld04_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp10_maxmiss100_ld04_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp10_maxmiss100_ld04_pve$pve)

# plot pca
prow_dp10_maxmiss100_ld04_pca_plot<-ggplot(prow_dp10_maxmiss100_ld04_pca, aes(PC1, PC2, col = prow_dp10_maxmiss100_ld04_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.4)) +
    xlab(paste0("PC1 (", signif(prow_dp10_maxmiss100_ld04_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp10_maxmiss100_ld04_pve$pve[2], 3), "%)"))

prow_dp10_maxmiss100_ld04_pca_plot





















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld05 ####
prow_dp10_maxmiss100_ld05_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_dp10_maxmiss100_ld05_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp10_maxmiss100_ld05_pca <- prow_dp10_maxmiss100_ld05_pca[,-1]
# set names
names(prow_dp10_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_dp10_maxmiss100_ld05_pca)[2:ncol(prow_dp10_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_dp10_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp10_maxmiss100_ld05_loc <- rep(NA, length(prow_dp10_maxmiss100_ld05_pca$ind))
prow_dp10_maxmiss100_ld05_loc[grep("VA", prow_dp10_maxmiss100_ld05_pca$ind)] <- "Virginia"
prow_dp10_maxmiss100_ld05_loc[grep("OH", prow_dp10_maxmiss100_ld05_pca$ind)] <- "Ohio"
prow_dp10_maxmiss100_ld05_loc[grep("AR", prow_dp10_maxmiss100_ld05_pca$ind)] <- "Arkansas"
prow_dp10_maxmiss100_ld05_loc[grep("WI", prow_dp10_maxmiss100_ld05_pca$ind)] <- "Wisconsin"
prow_dp10_maxmiss100_ld05_loc[grep("SC", prow_dp10_maxmiss100_ld05_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp10_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_dp10_maxmiss100_ld05_pca, prow_dp10_maxmiss100_ld05_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp10_maxmiss100_ld05_pve <- data.frame(PC = 1:20, pve = prow_dp10_maxmiss100_ld05_eigenval/sum(prow_dp10_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp10_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp10_maxmiss100_ld05_pve$pve)

# plot pca
prow_dp10_maxmiss100_ld05_pca_plot<-ggplot(prow_dp10_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_dp10_maxmiss100_ld05_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_dp10_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp10_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_dp10_maxmiss100_ld05_pca_plot















#### PROW DEPTH SET AT 15X ####

#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld0001 ####
prow_dp15_maxmiss100_ld0001_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld0001.eigenvec", col_names = FALSE)
prow_dp15_maxmiss100_ld0001_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld0001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp15_maxmiss100_ld0001_pca <- prow_dp15_maxmiss100_ld0001_pca[,-1]
# set names
names(prow_dp15_maxmiss100_ld0001_pca)[1] <- "ind"
names(prow_dp15_maxmiss100_ld0001_pca)[2:ncol(prow_dp15_maxmiss100_ld0001_pca)] <- paste0("PC", 1:(ncol(prow_dp15_maxmiss100_ld0001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp15_maxmiss100_ld0001_loc <- rep(NA, length(prow_dp15_maxmiss100_ld0001_pca$ind))
prow_dp15_maxmiss100_ld0001_loc[grep("VA", prow_dp15_maxmiss100_ld0001_pca$ind)] <- "Virginia"
prow_dp15_maxmiss100_ld0001_loc[grep("OH", prow_dp15_maxmiss100_ld0001_pca$ind)] <- "Ohio"
prow_dp15_maxmiss100_ld0001_loc[grep("AR", prow_dp15_maxmiss100_ld0001_pca$ind)] <- "Arkansas"
prow_dp15_maxmiss100_ld0001_loc[grep("WI", prow_dp15_maxmiss100_ld0001_pca$ind)] <- "Wisconsin"
prow_dp15_maxmiss100_ld0001_loc[grep("SC", prow_dp15_maxmiss100_ld0001_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp15_maxmiss100_ld0001_pca <- as_tibble(data.frame(prow_dp15_maxmiss100_ld0001_pca, prow_dp15_maxmiss100_ld0001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp15_maxmiss100_ld0001_pve <- data.frame(PC = 1:20, pve = prow_dp15_maxmiss100_ld0001_eigenval/sum(prow_dp15_maxmiss100_ld0001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp15_maxmiss100_ld0001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp15_maxmiss100_ld0001_pve$pve)

# plot pca
prow_dp15_maxmiss100_ld0001_pca_plot<-ggplot(prow_dp15_maxmiss100_ld0001_pca, aes(PC1, PC2, col = prow_dp15_maxmiss100_ld0001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.001)) +
    xlab(paste0("PC1 (", signif(prow_dp15_maxmiss100_ld0001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp15_maxmiss100_ld0001_pve$pve[2], 3), "%)"))

prow_dp15_maxmiss100_ld0001_pca_plot















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001 ####
prow_dp15_maxmiss100_ld001_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_dp15_maxmiss100_ld001_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp15_maxmiss100_ld001_pca <- prow_dp15_maxmiss100_ld001_pca[,-1]
# set names
names(prow_dp15_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_dp15_maxmiss100_ld001_pca)[2:ncol(prow_dp15_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_dp15_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp15_maxmiss100_ld001_loc <- rep(NA, length(prow_dp15_maxmiss100_ld001_pca$ind))
prow_dp15_maxmiss100_ld001_loc[grep("VA", prow_dp15_maxmiss100_ld001_pca$ind)] <- "Virginia"
prow_dp15_maxmiss100_ld001_loc[grep("OH", prow_dp15_maxmiss100_ld001_pca$ind)] <- "Ohio"
prow_dp15_maxmiss100_ld001_loc[grep("AR", prow_dp15_maxmiss100_ld001_pca$ind)] <- "Arkansas"
prow_dp15_maxmiss100_ld001_loc[grep("WI", prow_dp15_maxmiss100_ld001_pca$ind)] <- "Wisconsin"
prow_dp15_maxmiss100_ld001_loc[grep("SC", prow_dp15_maxmiss100_ld001_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp15_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_dp15_maxmiss100_ld001_pca, prow_dp15_maxmiss100_ld001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp15_maxmiss100_ld001_pve <- data.frame(PC = 1:20, pve = prow_dp15_maxmiss100_ld001_eigenval/sum(prow_dp15_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp15_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp15_maxmiss100_ld001_pve$pve)

# plot pca
prow_dp15_maxmiss100_ld001_pca_plot<-ggplot(prow_dp15_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_dp15_maxmiss100_ld001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_dp15_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp15_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_dp15_maxmiss100_ld001_pca_plot




















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01 ####
prow_dp15_maxmiss100_ld01_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.eigenvec", col_names = FALSE)
prow_dp15_maxmiss100_ld01_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld01.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp15_maxmiss100_ld01_pca <- prow_dp15_maxmiss100_ld01_pca[,-1]
# set names
names(prow_dp15_maxmiss100_ld01_pca)[1] <- "ind"
names(prow_dp15_maxmiss100_ld01_pca)[2:ncol(prow_dp15_maxmiss100_ld01_pca)] <- paste0("PC", 1:(ncol(prow_dp15_maxmiss100_ld01_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp15_maxmiss100_ld01_loc <- rep(NA, length(prow_dp15_maxmiss100_ld01_pca$ind))
prow_dp15_maxmiss100_ld01_loc[grep("VA", prow_dp15_maxmiss100_ld01_pca$ind)] <- "Virginia"
prow_dp15_maxmiss100_ld01_loc[grep("OH", prow_dp15_maxmiss100_ld01_pca$ind)] <- "Ohio"
prow_dp15_maxmiss100_ld01_loc[grep("AR", prow_dp15_maxmiss100_ld01_pca$ind)] <- "Arkansas"
prow_dp15_maxmiss100_ld01_loc[grep("WI", prow_dp15_maxmiss100_ld01_pca$ind)] <- "Wisconsin"
prow_dp15_maxmiss100_ld01_loc[grep("SC", prow_dp15_maxmiss100_ld01_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp15_maxmiss100_ld01_pca <- as_tibble(data.frame(prow_dp15_maxmiss100_ld01_pca, prow_dp15_maxmiss100_ld01_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp15_maxmiss100_ld01_pve <- data.frame(PC = 1:20, pve = prow_dp15_maxmiss100_ld01_eigenval/sum(prow_dp15_maxmiss100_ld01_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp15_maxmiss100_ld01_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp15_maxmiss100_ld01_pve$pve)

# plot pca
prow_dp15_maxmiss100_ld01_pca_plot<-ggplot(prow_dp15_maxmiss100_ld01_pca, aes(PC1, PC2, col = prow_dp15_maxmiss100_ld01_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.1)) +
    xlab(paste0("PC1 (", signif(prow_dp15_maxmiss100_ld01_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp15_maxmiss100_ld01_pve$pve[2], 3), "%)"))

prow_dp15_maxmiss100_ld01_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld02 ####
prow_dp15_maxmiss100_ld02_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld02.eigenvec", col_names = FALSE)
prow_dp15_maxmiss100_ld02_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld02.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp15_maxmiss100_ld02_pca <- prow_dp15_maxmiss100_ld02_pca[,-1]
# set names
names(prow_dp15_maxmiss100_ld02_pca)[1] <- "ind"
names(prow_dp15_maxmiss100_ld02_pca)[2:ncol(prow_dp15_maxmiss100_ld02_pca)] <- paste0("PC", 1:(ncol(prow_dp15_maxmiss100_ld02_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp15_maxmiss100_ld02_loc <- rep(NA, length(prow_dp15_maxmiss100_ld02_pca$ind))
prow_dp15_maxmiss100_ld02_loc[grep("VA", prow_dp15_maxmiss100_ld02_pca$ind)] <- "Virginia"
prow_dp15_maxmiss100_ld02_loc[grep("OH", prow_dp15_maxmiss100_ld02_pca$ind)] <- "Ohio"
prow_dp15_maxmiss100_ld02_loc[grep("AR", prow_dp15_maxmiss100_ld02_pca$ind)] <- "Arkansas"
prow_dp15_maxmiss100_ld02_loc[grep("WI", prow_dp15_maxmiss100_ld02_pca$ind)] <- "Wisconsin"
prow_dp15_maxmiss100_ld02_loc[grep("SC", prow_dp15_maxmiss100_ld02_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp15_maxmiss100_ld02_pca <- as_tibble(data.frame(prow_dp15_maxmiss100_ld02_pca, prow_dp15_maxmiss100_ld02_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp15_maxmiss100_ld02_pve <- data.frame(PC = 1:20, pve = prow_dp15_maxmiss100_ld02_eigenval/sum(prow_dp15_maxmiss100_ld02_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp15_maxmiss100_ld02_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp15_maxmiss100_ld02_pve$pve)

# plot pca
prow_dp15_maxmiss100_ld02_pca_plot<-ggplot(prow_dp15_maxmiss100_ld02_pca, aes(PC1, PC2, col = prow_dp15_maxmiss100_ld02_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.2)) +
    xlab(paste0("PC1 (", signif(prow_dp15_maxmiss100_ld02_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp15_maxmiss100_ld02_pve$pve[2], 3), "%)"))

prow_dp15_maxmiss100_ld02_pca_plot















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld03 ####
prow_dp15_maxmiss100_ld03_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld03.eigenvec", col_names = FALSE)
prow_dp15_maxmiss100_ld03_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld03.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp15_maxmiss100_ld03_pca <- prow_dp15_maxmiss100_ld03_pca[,-1]
# set names
names(prow_dp15_maxmiss100_ld03_pca)[1] <- "ind"
names(prow_dp15_maxmiss100_ld03_pca)[2:ncol(prow_dp15_maxmiss100_ld03_pca)] <- paste0("PC", 1:(ncol(prow_dp15_maxmiss100_ld03_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp15_maxmiss100_ld03_loc <- rep(NA, length(prow_dp15_maxmiss100_ld03_pca$ind))
prow_dp15_maxmiss100_ld03_loc[grep("VA", prow_dp15_maxmiss100_ld03_pca$ind)] <- "Virginia"
prow_dp15_maxmiss100_ld03_loc[grep("OH", prow_dp15_maxmiss100_ld03_pca$ind)] <- "Ohio"
prow_dp15_maxmiss100_ld03_loc[grep("AR", prow_dp15_maxmiss100_ld03_pca$ind)] <- "Arkansas"
prow_dp15_maxmiss100_ld03_loc[grep("WI", prow_dp15_maxmiss100_ld03_pca$ind)] <- "Wisconsin"
prow_dp15_maxmiss100_ld03_loc[grep("SC", prow_dp15_maxmiss100_ld03_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp15_maxmiss100_ld03_pca <- as_tibble(data.frame(prow_dp15_maxmiss100_ld03_pca, prow_dp15_maxmiss100_ld03_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp15_maxmiss100_ld03_pve <- data.frame(PC = 1:20, pve = prow_dp15_maxmiss100_ld03_eigenval/sum(prow_dp15_maxmiss100_ld03_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp15_maxmiss100_ld03_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp15_maxmiss100_ld03_pve$pve)

# plot pca
prow_dp15_maxmiss100_ld03_pca_plot<-ggplot(prow_dp15_maxmiss100_ld03_pca, aes(PC1, PC2, col = prow_dp15_maxmiss100_ld03_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.3)) +
    xlab(paste0("PC1 (", signif(prow_dp15_maxmiss100_ld03_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp15_maxmiss100_ld03_pve$pve[2], 3), "%)"))

prow_dp15_maxmiss100_ld03_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld04 ####
prow_dp15_maxmiss100_ld04_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld04.eigenvec", col_names = FALSE)
prow_dp15_maxmiss100_ld04_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld04.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp15_maxmiss100_ld04_pca <- prow_dp15_maxmiss100_ld04_pca[,-1]
# set names
names(prow_dp15_maxmiss100_ld04_pca)[1] <- "ind"
names(prow_dp15_maxmiss100_ld04_pca)[2:ncol(prow_dp15_maxmiss100_ld04_pca)] <- paste0("PC", 1:(ncol(prow_dp15_maxmiss100_ld04_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp15_maxmiss100_ld04_loc <- rep(NA, length(prow_dp15_maxmiss100_ld04_pca$ind))
prow_dp15_maxmiss100_ld04_loc[grep("VA", prow_dp15_maxmiss100_ld04_pca$ind)] <- "Virginia"
prow_dp15_maxmiss100_ld04_loc[grep("OH", prow_dp15_maxmiss100_ld04_pca$ind)] <- "Ohio"
prow_dp15_maxmiss100_ld04_loc[grep("AR", prow_dp15_maxmiss100_ld04_pca$ind)] <- "Arkansas"
prow_dp15_maxmiss100_ld04_loc[grep("WI", prow_dp15_maxmiss100_ld04_pca$ind)] <- "Wisconsin"
prow_dp15_maxmiss100_ld04_loc[grep("SC", prow_dp15_maxmiss100_ld04_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp15_maxmiss100_ld04_pca <- as_tibble(data.frame(prow_dp15_maxmiss100_ld04_pca, prow_dp15_maxmiss100_ld04_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp15_maxmiss100_ld04_pve <- data.frame(PC = 1:20, pve = prow_dp15_maxmiss100_ld04_eigenval/sum(prow_dp15_maxmiss100_ld04_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp15_maxmiss100_ld04_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp15_maxmiss100_ld04_pve$pve)

# plot pca
prow_dp15_maxmiss100_ld04_pca_plot<-ggplot(prow_dp15_maxmiss100_ld04_pca, aes(PC1, PC2, col = prow_dp15_maxmiss100_ld04_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.4)) +
    xlab(paste0("PC1 (", signif(prow_dp15_maxmiss100_ld04_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp15_maxmiss100_ld04_pve$pve[2], 3), "%)"))

prow_dp15_maxmiss100_ld04_pca_plot





















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld05 ####
prow_dp15_maxmiss100_ld05_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_dp15_maxmiss100_ld05_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp15_maxmiss100_ld05_pca <- prow_dp15_maxmiss100_ld05_pca[,-1]
# set names
names(prow_dp15_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_dp15_maxmiss100_ld05_pca)[2:ncol(prow_dp15_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_dp15_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp15_maxmiss100_ld05_loc <- rep(NA, length(prow_dp15_maxmiss100_ld05_pca$ind))
prow_dp15_maxmiss100_ld05_loc[grep("VA", prow_dp15_maxmiss100_ld05_pca$ind)] <- "Virginia"
prow_dp15_maxmiss100_ld05_loc[grep("OH", prow_dp15_maxmiss100_ld05_pca$ind)] <- "Ohio"
prow_dp15_maxmiss100_ld05_loc[grep("AR", prow_dp15_maxmiss100_ld05_pca$ind)] <- "Arkansas"
prow_dp15_maxmiss100_ld05_loc[grep("WI", prow_dp15_maxmiss100_ld05_pca$ind)] <- "Wisconsin"
prow_dp15_maxmiss100_ld05_loc[grep("SC", prow_dp15_maxmiss100_ld05_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp15_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_dp15_maxmiss100_ld05_pca, prow_dp15_maxmiss100_ld05_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp15_maxmiss100_ld05_pve <- data.frame(PC = 1:20, pve = prow_dp15_maxmiss100_ld05_eigenval/sum(prow_dp15_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp15_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp15_maxmiss100_ld05_pve$pve)

# plot pca
prow_dp15_maxmiss100_ld05_pca_plot<-ggplot(prow_dp15_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_dp15_maxmiss100_ld05_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_dp15_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp15_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_dp15_maxmiss100_ld05_pca_plot







#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld08 ####
prow_dp15_maxmiss100_ld08_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld08.eigenvec", col_names = FALSE)
prow_dp15_maxmiss100_ld08_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp15_maxmiss100_ld08.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp15_maxmiss100_ld08_pca <- prow_dp15_maxmiss100_ld08_pca[,-1]
# set names
names(prow_dp15_maxmiss100_ld08_pca)[1] <- "ind"
names(prow_dp15_maxmiss100_ld08_pca)[2:ncol(prow_dp15_maxmiss100_ld08_pca)] <- paste0("PC", 1:(ncol(prow_dp15_maxmiss100_ld08_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp15_maxmiss100_ld08_loc <- rep(NA, length(prow_dp15_maxmiss100_ld08_pca$ind))
prow_dp15_maxmiss100_ld08_loc[grep("VA", prow_dp15_maxmiss100_ld08_pca$ind)] <- "Virginia"
prow_dp15_maxmiss100_ld08_loc[grep("OH", prow_dp15_maxmiss100_ld08_pca$ind)] <- "Ohio"
prow_dp15_maxmiss100_ld08_loc[grep("AR", prow_dp15_maxmiss100_ld08_pca$ind)] <- "Arkansas"
prow_dp15_maxmiss100_ld08_loc[grep("WI", prow_dp15_maxmiss100_ld08_pca$ind)] <- "Wisconsin"
prow_dp15_maxmiss100_ld08_loc[grep("SC", prow_dp15_maxmiss100_ld08_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp15_maxmiss100_ld08_pca <- as_tibble(data.frame(prow_dp15_maxmiss100_ld08_pca, prow_dp15_maxmiss100_ld08_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp15_maxmiss100_ld08_pve <- data.frame(PC = 1:20, pve = prow_dp15_maxmiss100_ld08_eigenval/sum(prow_dp15_maxmiss100_ld08_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp15_maxmiss100_ld08_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp15_maxmiss100_ld08_pve$pve)

# plot pca
prow_dp15_maxmiss100_ld08_pca_plot<-ggplot(prow_dp15_maxmiss100_ld08_pca, aes(PC1, PC2, col = prow_dp15_maxmiss100_ld08_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.8)) +
    xlab(paste0("PC1 (", signif(prow_dp15_maxmiss100_ld08_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp15_maxmiss100_ld08_pve$pve[2], 3), "%)"))

prow_dp15_maxmiss100_ld08_pca_plot

















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld08 ####
prow_dp10_maxmiss100_ld08_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld08.eigenvec", col_names = FALSE)
prow_dp10_maxmiss100_ld08_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp10_maxmiss100_ld08.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp10_maxmiss100_ld08_pca <- prow_dp10_maxmiss100_ld08_pca[,-1]
# set names
names(prow_dp10_maxmiss100_ld08_pca)[1] <- "ind"
names(prow_dp10_maxmiss100_ld08_pca)[2:ncol(prow_dp10_maxmiss100_ld08_pca)] <- paste0("PC", 1:(ncol(prow_dp10_maxmiss100_ld08_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp10_maxmiss100_ld08_loc <- rep(NA, length(prow_dp10_maxmiss100_ld08_pca$ind))
prow_dp10_maxmiss100_ld08_loc[grep("VA", prow_dp10_maxmiss100_ld08_pca$ind)] <- "Virginia"
prow_dp10_maxmiss100_ld08_loc[grep("OH", prow_dp10_maxmiss100_ld08_pca$ind)] <- "Ohio"
prow_dp10_maxmiss100_ld08_loc[grep("AR", prow_dp10_maxmiss100_ld08_pca$ind)] <- "Arkansas"
prow_dp10_maxmiss100_ld08_loc[grep("WI", prow_dp10_maxmiss100_ld08_pca$ind)] <- "Wisconsin"
prow_dp10_maxmiss100_ld08_loc[grep("SC", prow_dp10_maxmiss100_ld08_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp10_maxmiss100_ld08_pca <- as_tibble(data.frame(prow_dp10_maxmiss100_ld08_pca, prow_dp10_maxmiss100_ld08_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp10_maxmiss100_ld08_pve <- data.frame(PC = 1:20, pve = prow_dp10_maxmiss100_ld08_eigenval/sum(prow_dp10_maxmiss100_ld08_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp10_maxmiss100_ld08_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp10_maxmiss100_ld08_pve$pve)

# plot pca
prow_dp10_maxmiss100_ld08_pca_plot<-ggplot(prow_dp10_maxmiss100_ld08_pca, aes(PC1, PC2, col = prow_dp10_maxmiss100_ld08_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.8)) +
    xlab(paste0("PC1 (", signif(prow_dp10_maxmiss100_ld08_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp10_maxmiss100_ld08_pve$pve[2], 3), "%)"))

prow_dp10_maxmiss100_ld08_pca_plot


















#### PROW DEPTH SET AT 20X ####

#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld0001 ####
prow_dp20_maxmiss100_ld0001_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld0001.eigenvec", col_names = FALSE)
prow_dp20_maxmiss100_ld0001_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld0001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp20_maxmiss100_ld0001_pca <- prow_dp20_maxmiss100_ld0001_pca[,-1]
# set names
names(prow_dp20_maxmiss100_ld0001_pca)[1] <- "ind"
names(prow_dp20_maxmiss100_ld0001_pca)[2:ncol(prow_dp20_maxmiss100_ld0001_pca)] <- paste0("PC", 1:(ncol(prow_dp20_maxmiss100_ld0001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp20_maxmiss100_ld0001_loc <- rep(NA, length(prow_dp20_maxmiss100_ld0001_pca$ind))
prow_dp20_maxmiss100_ld0001_loc[grep("VA", prow_dp20_maxmiss100_ld0001_pca$ind)] <- "Virginia"
prow_dp20_maxmiss100_ld0001_loc[grep("OH", prow_dp20_maxmiss100_ld0001_pca$ind)] <- "Ohio"
prow_dp20_maxmiss100_ld0001_loc[grep("AR", prow_dp20_maxmiss100_ld0001_pca$ind)] <- "Arkansas"
prow_dp20_maxmiss100_ld0001_loc[grep("WI", prow_dp20_maxmiss100_ld0001_pca$ind)] <- "Wisconsin"
prow_dp20_maxmiss100_ld0001_loc[grep("SC", prow_dp20_maxmiss100_ld0001_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp20_maxmiss100_ld0001_pca <- as_tibble(data.frame(prow_dp20_maxmiss100_ld0001_pca, prow_dp20_maxmiss100_ld0001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp20_maxmiss100_ld0001_pve <- data.frame(PC = 1:20, pve = prow_dp20_maxmiss100_ld0001_eigenval/sum(prow_dp20_maxmiss100_ld0001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp20_maxmiss100_ld0001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp20_maxmiss100_ld0001_pve$pve)

# plot pca
prow_dp20_maxmiss100_ld0001_pca_plot<-ggplot(prow_dp20_maxmiss100_ld0001_pca, aes(PC1, PC2, col = prow_dp20_maxmiss100_ld0001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.001)) +
    xlab(paste0("PC1 (", signif(prow_dp20_maxmiss100_ld0001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp20_maxmiss100_ld0001_pve$pve[2], 3), "%)"))

prow_dp20_maxmiss100_ld0001_pca_plot















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld001 ####
prow_dp20_maxmiss100_ld001_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_dp20_maxmiss100_ld001_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp20_maxmiss100_ld001_pca <- prow_dp20_maxmiss100_ld001_pca[,-1]
# set names
names(prow_dp20_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_dp20_maxmiss100_ld001_pca)[2:ncol(prow_dp20_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_dp20_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp20_maxmiss100_ld001_loc <- rep(NA, length(prow_dp20_maxmiss100_ld001_pca$ind))
prow_dp20_maxmiss100_ld001_loc[grep("VA", prow_dp20_maxmiss100_ld001_pca$ind)] <- "Virginia"
prow_dp20_maxmiss100_ld001_loc[grep("OH", prow_dp20_maxmiss100_ld001_pca$ind)] <- "Ohio"
prow_dp20_maxmiss100_ld001_loc[grep("AR", prow_dp20_maxmiss100_ld001_pca$ind)] <- "Arkansas"
prow_dp20_maxmiss100_ld001_loc[grep("WI", prow_dp20_maxmiss100_ld001_pca$ind)] <- "Wisconsin"
prow_dp20_maxmiss100_ld001_loc[grep("SC", prow_dp20_maxmiss100_ld001_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp20_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_dp20_maxmiss100_ld001_pca, prow_dp20_maxmiss100_ld001_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp20_maxmiss100_ld001_pve <- data.frame(PC = 1:20, pve = prow_dp20_maxmiss100_ld001_eigenval/sum(prow_dp20_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp20_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp20_maxmiss100_ld001_pve$pve)

# plot pca
prow_dp20_maxmiss100_ld001_pca_plot<-ggplot(prow_dp20_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_dp20_maxmiss100_ld001_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_dp20_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp20_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_dp20_maxmiss100_ld001_pca_plot




















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld01 ####
prow_dp20_maxmiss100_ld01_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld01.eigenvec", col_names = FALSE)
prow_dp20_maxmiss100_ld01_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld01.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp20_maxmiss100_ld01_pca <- prow_dp20_maxmiss100_ld01_pca[,-1]
# set names
names(prow_dp20_maxmiss100_ld01_pca)[1] <- "ind"
names(prow_dp20_maxmiss100_ld01_pca)[2:ncol(prow_dp20_maxmiss100_ld01_pca)] <- paste0("PC", 1:(ncol(prow_dp20_maxmiss100_ld01_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp20_maxmiss100_ld01_loc <- rep(NA, length(prow_dp20_maxmiss100_ld01_pca$ind))
prow_dp20_maxmiss100_ld01_loc[grep("VA", prow_dp20_maxmiss100_ld01_pca$ind)] <- "Virginia"
prow_dp20_maxmiss100_ld01_loc[grep("OH", prow_dp20_maxmiss100_ld01_pca$ind)] <- "Ohio"
prow_dp20_maxmiss100_ld01_loc[grep("AR", prow_dp20_maxmiss100_ld01_pca$ind)] <- "Arkansas"
prow_dp20_maxmiss100_ld01_loc[grep("WI", prow_dp20_maxmiss100_ld01_pca$ind)] <- "Wisconsin"
prow_dp20_maxmiss100_ld01_loc[grep("SC", prow_dp20_maxmiss100_ld01_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp20_maxmiss100_ld01_pca <- as_tibble(data.frame(prow_dp20_maxmiss100_ld01_pca, prow_dp20_maxmiss100_ld01_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp20_maxmiss100_ld01_pve <- data.frame(PC = 1:20, pve = prow_dp20_maxmiss100_ld01_eigenval/sum(prow_dp20_maxmiss100_ld01_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp20_maxmiss100_ld01_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp20_maxmiss100_ld01_pve$pve)

# plot pca
prow_dp20_maxmiss100_ld01_pca_plot<-ggplot(prow_dp20_maxmiss100_ld01_pca, aes(PC1, PC2, col = prow_dp20_maxmiss100_ld01_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.1)) +
    xlab(paste0("PC1 (", signif(prow_dp20_maxmiss100_ld01_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp20_maxmiss100_ld01_pve$pve[2], 3), "%)"))

prow_dp20_maxmiss100_ld01_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld02 ####
prow_dp20_maxmiss100_ld02_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld02.eigenvec", col_names = FALSE)
prow_dp20_maxmiss100_ld02_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld02.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp20_maxmiss100_ld02_pca <- prow_dp20_maxmiss100_ld02_pca[,-1]
# set names
names(prow_dp20_maxmiss100_ld02_pca)[1] <- "ind"
names(prow_dp20_maxmiss100_ld02_pca)[2:ncol(prow_dp20_maxmiss100_ld02_pca)] <- paste0("PC", 1:(ncol(prow_dp20_maxmiss100_ld02_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp20_maxmiss100_ld02_loc <- rep(NA, length(prow_dp20_maxmiss100_ld02_pca$ind))
prow_dp20_maxmiss100_ld02_loc[grep("VA", prow_dp20_maxmiss100_ld02_pca$ind)] <- "Virginia"
prow_dp20_maxmiss100_ld02_loc[grep("OH", prow_dp20_maxmiss100_ld02_pca$ind)] <- "Ohio"
prow_dp20_maxmiss100_ld02_loc[grep("AR", prow_dp20_maxmiss100_ld02_pca$ind)] <- "Arkansas"
prow_dp20_maxmiss100_ld02_loc[grep("WI", prow_dp20_maxmiss100_ld02_pca$ind)] <- "Wisconsin"
prow_dp20_maxmiss100_ld02_loc[grep("SC", prow_dp20_maxmiss100_ld02_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp20_maxmiss100_ld02_pca <- as_tibble(data.frame(prow_dp20_maxmiss100_ld02_pca, prow_dp20_maxmiss100_ld02_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp20_maxmiss100_ld02_pve <- data.frame(PC = 1:20, pve = prow_dp20_maxmiss100_ld02_eigenval/sum(prow_dp20_maxmiss100_ld02_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp20_maxmiss100_ld02_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp20_maxmiss100_ld02_pve$pve)

# plot pca
prow_dp20_maxmiss100_ld02_pca_plot<-ggplot(prow_dp20_maxmiss100_ld02_pca, aes(PC1, PC2, col = prow_dp20_maxmiss100_ld02_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.2)) +
    xlab(paste0("PC1 (", signif(prow_dp20_maxmiss100_ld02_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp20_maxmiss100_ld02_pve$pve[2], 3), "%)"))

prow_dp20_maxmiss100_ld02_pca_plot















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld03 ####
prow_dp20_maxmiss100_ld03_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld03.eigenvec", col_names = FALSE)
prow_dp20_maxmiss100_ld03_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld03.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp20_maxmiss100_ld03_pca <- prow_dp20_maxmiss100_ld03_pca[,-1]
# set names
names(prow_dp20_maxmiss100_ld03_pca)[1] <- "ind"
names(prow_dp20_maxmiss100_ld03_pca)[2:ncol(prow_dp20_maxmiss100_ld03_pca)] <- paste0("PC", 1:(ncol(prow_dp20_maxmiss100_ld03_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp20_maxmiss100_ld03_loc <- rep(NA, length(prow_dp20_maxmiss100_ld03_pca$ind))
prow_dp20_maxmiss100_ld03_loc[grep("VA", prow_dp20_maxmiss100_ld03_pca$ind)] <- "Virginia"
prow_dp20_maxmiss100_ld03_loc[grep("OH", prow_dp20_maxmiss100_ld03_pca$ind)] <- "Ohio"
prow_dp20_maxmiss100_ld03_loc[grep("AR", prow_dp20_maxmiss100_ld03_pca$ind)] <- "Arkansas"
prow_dp20_maxmiss100_ld03_loc[grep("WI", prow_dp20_maxmiss100_ld03_pca$ind)] <- "Wisconsin"
prow_dp20_maxmiss100_ld03_loc[grep("SC", prow_dp20_maxmiss100_ld03_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp20_maxmiss100_ld03_pca <- as_tibble(data.frame(prow_dp20_maxmiss100_ld03_pca, prow_dp20_maxmiss100_ld03_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp20_maxmiss100_ld03_pve <- data.frame(PC = 1:20, pve = prow_dp20_maxmiss100_ld03_eigenval/sum(prow_dp20_maxmiss100_ld03_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp20_maxmiss100_ld03_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp20_maxmiss100_ld03_pve$pve)

# plot pca
prow_dp20_maxmiss100_ld03_pca_plot<-ggplot(prow_dp20_maxmiss100_ld03_pca, aes(PC1, PC2, col = prow_dp20_maxmiss100_ld03_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.3)) +
    xlab(paste0("PC1 (", signif(prow_dp20_maxmiss100_ld03_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp20_maxmiss100_ld03_pve$pve[2], 3), "%)"))

prow_dp20_maxmiss100_ld03_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld04 ####
prow_dp20_maxmiss100_ld04_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld04.eigenvec", col_names = FALSE)
prow_dp20_maxmiss100_ld04_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld04.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp20_maxmiss100_ld04_pca <- prow_dp20_maxmiss100_ld04_pca[,-1]
# set names
names(prow_dp20_maxmiss100_ld04_pca)[1] <- "ind"
names(prow_dp20_maxmiss100_ld04_pca)[2:ncol(prow_dp20_maxmiss100_ld04_pca)] <- paste0("PC", 1:(ncol(prow_dp20_maxmiss100_ld04_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp20_maxmiss100_ld04_loc <- rep(NA, length(prow_dp20_maxmiss100_ld04_pca$ind))
prow_dp20_maxmiss100_ld04_loc[grep("VA", prow_dp20_maxmiss100_ld04_pca$ind)] <- "Virginia"
prow_dp20_maxmiss100_ld04_loc[grep("OH", prow_dp20_maxmiss100_ld04_pca$ind)] <- "Ohio"
prow_dp20_maxmiss100_ld04_loc[grep("AR", prow_dp20_maxmiss100_ld04_pca$ind)] <- "Arkansas"
prow_dp20_maxmiss100_ld04_loc[grep("WI", prow_dp20_maxmiss100_ld04_pca$ind)] <- "Wisconsin"
prow_dp20_maxmiss100_ld04_loc[grep("SC", prow_dp20_maxmiss100_ld04_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp20_maxmiss100_ld04_pca <- as_tibble(data.frame(prow_dp20_maxmiss100_ld04_pca, prow_dp20_maxmiss100_ld04_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp20_maxmiss100_ld04_pve <- data.frame(PC = 1:20, pve = prow_dp20_maxmiss100_ld04_eigenval/sum(prow_dp20_maxmiss100_ld04_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp20_maxmiss100_ld04_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp20_maxmiss100_ld04_pve$pve)

# plot pca
prow_dp20_maxmiss100_ld04_pca_plot<-ggplot(prow_dp20_maxmiss100_ld04_pca, aes(PC1, PC2, col = prow_dp20_maxmiss100_ld04_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.4)) +
    xlab(paste0("PC1 (", signif(prow_dp20_maxmiss100_ld04_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp20_maxmiss100_ld04_pve$pve[2], 3), "%)"))

prow_dp20_maxmiss100_ld04_pca_plot





















#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld05 ####
prow_dp20_maxmiss100_ld05_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_dp20_maxmiss100_ld05_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp20_maxmiss100_ld05_pca <- prow_dp20_maxmiss100_ld05_pca[,-1]
# set names
names(prow_dp20_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_dp20_maxmiss100_ld05_pca)[2:ncol(prow_dp20_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_dp20_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp20_maxmiss100_ld05_loc <- rep(NA, length(prow_dp20_maxmiss100_ld05_pca$ind))
prow_dp20_maxmiss100_ld05_loc[grep("VA", prow_dp20_maxmiss100_ld05_pca$ind)] <- "Virginia"
prow_dp20_maxmiss100_ld05_loc[grep("OH", prow_dp20_maxmiss100_ld05_pca$ind)] <- "Ohio"
prow_dp20_maxmiss100_ld05_loc[grep("AR", prow_dp20_maxmiss100_ld05_pca$ind)] <- "Arkansas"
prow_dp20_maxmiss100_ld05_loc[grep("WI", prow_dp20_maxmiss100_ld05_pca$ind)] <- "Wisconsin"
prow_dp20_maxmiss100_ld05_loc[grep("SC", prow_dp20_maxmiss100_ld05_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp20_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_dp20_maxmiss100_ld05_pca, prow_dp20_maxmiss100_ld05_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp20_maxmiss100_ld05_pve <- data.frame(PC = 1:20, pve = prow_dp20_maxmiss100_ld05_eigenval/sum(prow_dp20_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp20_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp20_maxmiss100_ld05_pve$pve)

# plot pca
prow_dp20_maxmiss100_ld05_pca_plot<-ggplot(prow_dp20_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_dp20_maxmiss100_ld05_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_dp20_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp20_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_dp20_maxmiss100_ld05_pca_plot












#### ++ PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld08 ####
prow_dp20_maxmiss100_ld08_pca <- read_table("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld08.eigenvec", col_names = FALSE)
prow_dp20_maxmiss100_ld08_eigenval <- scan("./inputs/PROW_FILTERED_renamed_q30_minac1_maf05_dp20_maxmiss100_ld08.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_dp20_maxmiss100_ld08_pca <- prow_dp20_maxmiss100_ld08_pca[,-1]
# set names
names(prow_dp20_maxmiss100_ld08_pca)[1] <- "ind"
names(prow_dp20_maxmiss100_ld08_pca)[2:ncol(prow_dp20_maxmiss100_ld08_pca)] <- paste0("PC", 1:(ncol(prow_dp20_maxmiss100_ld08_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops
# location
prow_dp20_maxmiss100_ld08_loc <- rep(NA, length(prow_dp20_maxmiss100_ld08_pca$ind))
prow_dp20_maxmiss100_ld08_loc[grep("VA", prow_dp20_maxmiss100_ld08_pca$ind)] <- "Virginia"
prow_dp20_maxmiss100_ld08_loc[grep("OH", prow_dp20_maxmiss100_ld08_pca$ind)] <- "Ohio"
prow_dp20_maxmiss100_ld08_loc[grep("AR", prow_dp20_maxmiss100_ld08_pca$ind)] <- "Arkansas"
prow_dp20_maxmiss100_ld08_loc[grep("WI", prow_dp20_maxmiss100_ld08_pca$ind)] <- "Wisconsin"
prow_dp20_maxmiss100_ld08_loc[grep("SC", prow_dp20_maxmiss100_ld08_pca$ind)] <- "South Carolina"

# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_dp20_maxmiss100_ld08_pca <- as_tibble(data.frame(prow_dp20_maxmiss100_ld08_pca, prow_dp20_maxmiss100_ld08_loc))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_dp20_maxmiss100_ld08_pve <- data.frame(PC = 1:20, pve = prow_dp20_maxmiss100_ld08_eigenval/sum(prow_dp20_maxmiss100_ld08_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_dp20_maxmiss100_ld08_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_dp20_maxmiss100_ld08_pve$pve)

# plot pca
prow_dp20_maxmiss100_ld08_pca_plot<-ggplot(prow_dp20_maxmiss100_ld08_pca, aes(PC1, PC2, col = prow_dp20_maxmiss100_ld08_loc)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("darkcyan", "goldenrod2", "darkslategrey", "sienna2", "lightblue3")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Location (A. protonotaria)") +
    theme(legend.title = element_text(face = "italic")) +
    ggtitle(bquote('LD'~r^2:0.8)) +
    xlab(paste0("PC1 (", signif(prow_dp20_maxmiss100_ld08_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_dp20_maxmiss100_ld08_pve$pve[2], 3), "%)"))

prow_dp20_maxmiss100_ld08_pca_plot


















#### ++++ ALL PLOTS 5X, variable LD thresholds ####
ggarrange(
    cerw_dp05_maxmiss100_ld0001_pca_plot, cerw_dp05_maxmiss100_ld001_pca_plot, cerw_dp05_maxmiss100_ld01_pca_plot, cerw_dp05_maxmiss100_ld02_pca_plot, cerw_dp05_maxmiss100_ld03_pca_plot, cerw_dp05_maxmiss100_ld04_pca_plot, cerw_dp05_maxmiss100_ld05_pca_plot, cerw_dp05_maxmiss100_ld08_pca_plot,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    ncol = 4, nrow = 2, common.legend = TRUE)

# 8x16, landscape


ggarrange(
    prow_dp05_maxmiss100_ld0001_pca_plot, prow_dp05_maxmiss100_ld001_pca_plot, prow_dp05_maxmiss100_ld01_pca_plot, prow_dp05_maxmiss100_ld02_pca_plot, prow_dp05_maxmiss100_ld03_pca_plot, prow_dp05_maxmiss100_ld04_pca_plot, prow_dp05_maxmiss100_ld05_pca_plot, prow_dp05_maxmiss100_ld08_pca_plot,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    ncol = 4, nrow = 2, common.legend = TRUE)

# 8x16, landscape








#### ++++ ALL PLOTS 10X, variable LD thresholds ####
ggarrange(
    cerw_dp10_maxmiss100_ld0001_pca_plot, cerw_dp10_maxmiss100_ld001_pca_plot, cerw_dp10_maxmiss100_ld01_pca_plot, cerw_dp10_maxmiss100_ld02_pca_plot, cerw_dp10_maxmiss100_ld03_pca_plot, cerw_dp10_maxmiss100_ld04_pca_plot, cerw_dp10_maxmiss100_ld05_pca_plot, cerw_dp10_maxmiss100_ld08_pca_plot,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    ncol = 4, nrow = 2, common.legend = TRUE)

# 8x16, landscape


ggarrange(
    prow_dp10_maxmiss100_ld0001_pca_plot, prow_dp10_maxmiss100_ld001_pca_plot, prow_dp10_maxmiss100_ld01_pca_plot, prow_dp10_maxmiss100_ld02_pca_plot, prow_dp10_maxmiss100_ld03_pca_plot, prow_dp10_maxmiss100_ld04_pca_plot, prow_dp10_maxmiss100_ld05_pca_plot, prow_dp10_maxmiss100_ld08_pca_plot,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    ncol = 4, nrow = 2, common.legend = TRUE)

# 8x16, landscape












#### ++++ ALL PLOTS 15X, variable LD thresholds ####
ggarrange(
    cerw_dp15_maxmiss100_ld0001_pca_plot, cerw_dp15_maxmiss100_ld001_pca_plot, cerw_dp15_maxmiss100_ld01_pca_plot, cerw_dp15_maxmiss100_ld02_pca_plot, cerw_dp15_maxmiss100_ld03_pca_plot, cerw_dp15_maxmiss100_ld04_pca_plot, cerw_dp15_maxmiss100_ld05_pca_plot, cerw_dp15_maxmiss100_ld08_pca_plot,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    ncol = 4, nrow = 2, common.legend = TRUE)

# 8x16, landscape


ggarrange(
    prow_dp15_maxmiss100_ld0001_pca_plot, prow_dp15_maxmiss100_ld001_pca_plot, prow_dp15_maxmiss100_ld01_pca_plot, prow_dp15_maxmiss100_ld02_pca_plot, prow_dp15_maxmiss100_ld03_pca_plot, prow_dp15_maxmiss100_ld04_pca_plot, prow_dp15_maxmiss100_ld05_pca_plot, prow_dp15_maxmiss100_ld08_pca_plot,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    ncol = 4, nrow = 2, common.legend = TRUE)

# 8x16, landscape







#### ++++ ALL PLOTS 20X, variable LD thresholds ####
ggarrange(
    cerw_dp20_maxmiss100_ld0001_pca_plot, cerw_dp20_maxmiss100_ld001_pca_plot, cerw_dp20_maxmiss100_ld01_pca_plot, cerw_dp20_maxmiss100_ld02_pca_plot, cerw_dp20_maxmiss100_ld03_pca_plot, cerw_dp20_maxmiss100_ld04_pca_plot, cerw_dp20_maxmiss100_ld05_pca_plot, cerw_dp20_maxmiss100_ld08_pca_plot,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    ncol = 4, nrow = 2, common.legend = TRUE)

# 8x16, landscape


ggarrange(
    prow_dp20_maxmiss100_ld0001_pca_plot, prow_dp20_maxmiss100_ld001_pca_plot, prow_dp20_maxmiss100_ld01_pca_plot, prow_dp20_maxmiss100_ld02_pca_plot, prow_dp20_maxmiss100_ld03_pca_plot, prow_dp20_maxmiss100_ld04_pca_plot, prow_dp20_maxmiss100_ld05_pca_plot, prow_dp20_maxmiss100_ld08_pca_plot,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    ncol = 4, nrow = 2, common.legend = TRUE)

# 8x16, landscape
