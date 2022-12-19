#### Number of SNPs versus LD threshold ####

setwd("~/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_snp/03_SNP_20220513")

#### LOAD LIBRARIES ####

library(ggplot2)
library(scales)

data <- read.csv("NumSNPs_LD.csv", header=TRUE)
data$LD_threshold <- as.factor(data$LD_threshold)
str(data)

data_10x<-subset(data, Depth =="10X")
levels(data_10x$Depth)
str(data_10x)
data_10x$Depth <- factor(data_10x$Depth) # factor it again to remove 15X
levels(data_10x$Depth)


data_15x<-subset(data, Depth =="15X")
levels(data_15x$Depth)
str(data_15x)
data_15x$Depth <- factor(data_15x$Depth) # factor it again to remove 15X
levels(data_15x$Depth)


ggplot(data_10x, aes(x = LD_threshold, y = Number_of_SNPs)) +
    geom_point(aes(col=Species), size = 3) +
    labs(color = "Host Species") +
    ylab("Number of SNPs") +
    xlab("LD threshold") +
    ggtitle("DP = 10X") +
    scale_y_continuous(labels = scales::comma) +
    theme_bw()


# A6, landscape



ggplot(data_15x, aes(x = LD_threshold, y = Number_of_SNPs)) +
    geom_point(aes(col=Species), size = 3) +
    labs(color = "Host Species") +
    ylab("Number of SNPs") +
    xlab("LD threshold") +
    ggtitle("DP = 15X") +
    scale_y_continuous(labels = scales::comma) +
    theme_bw()

#  A6, landscape

