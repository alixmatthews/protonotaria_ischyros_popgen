## PLINK analyses with the 3 PROW samples removed

### Files: 

`04_plink_LD_PCA_$SPP_reducedref.slurm` - this is with 5X and 10X with 50%, 75%, and 100% missingness thresholds (LD r^2 value is 0.5)

`04_plink_LD_PCA_$SPP_1520x_reducedref.slurm` - this is with 15X and 20X with 50%, 75%, and 100% missingness thresholds. Decided to do this one also, but just at a separate time (LD r^2 value is 0.5)

`04_plink_LD08_PCA_$SPP_reducedref.slurm` - this is with 5X, 10X, 15X, and 20X with 50%, 75%, and 100% missingness thresholds (LD r^2 value is 0.8)

`04_plink_LD02_PCA_$SPP_reducedref.slurm` - this is with 5X, 10X, 15X, and 20X with 50%, 75%, and 100% missingness thresholds (LD r^2 value is 0.2)

`plink_PCA.R` - R script to make PCAs from files in the `inputs` directory (.eigenval and .eigenvec files). Outputs from this R script are in the `output` directory
