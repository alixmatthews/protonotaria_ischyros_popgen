## Let's try an IBD analysis on these data

Isolation-by-distance analysis of the pairwise Fst per sample and the pairwise geographic distance between sampling sites

---

**1.** Geographic distance matrices

**1.1.** Run Geographic Distance Matrix Generator

- installed Geographic distance matrix generator (https://biodiversityinformatics.amnh.org/open_source/gdmg/)
- extracted files
- moved files into /home/lse305
- Run GDMG 
```
cd /home/lse305/GeographicDistanceMatrixGenerator_v1.2.3/
./Start\ GDMG\ Mac_Linux.sh 
```

Inputs for GDMG:
- each sample has its own GPS location, I just slightly changed the GPS points (e.g., 1/1000th degree different) for samples that didn't have unique GPS coordinates)
  - `CERW_GDM_input.txt`
  - `PROW_GDM_input.txt`

Outputs created by GDMG:
  - `CERW_GDM_output.txt`
  - `PROW_GDM_output.txt`

---

**2.** Run IBD analyses in R

**2.1.** R script
- This has some other instructions within (i.e., where the Fst input comes from) 
- This also runs Mantel tests to see the correlation between Fst matrices and geographic distance matrices
- Also creates figures

Figures:
  - `CERW_IBD.pdf`
  - `PROW_IBD.pdf`
  - `CERW_PROW_IBD.pdf` - results together



