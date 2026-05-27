# Principal Component and Ordination Analyses

**Input dataset**: Unlinked SNP dataset in PLINK format (`biallelic_snps_thinned.bed/.bim/.fam`)

---

## Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[N_SNPS]` | SNP count in final analysis VCF |
| `[N_SAMPLES]` | Number of individuals |
| `[PCAONE_N_PC]` | `config["parameters"]["PCAone"]["PCnum"]` (default: 10) |
| `[PCAONE_SVD]` | `config["parameters"]["PCAone"]["SVD_method"]` (default: 3 = full SVD) |
| `[EMU_N_PC]` | `config["parameters"]["PCAone"]["PCnum"]` (same parameter) |
| `[PCOA_DISTANCE]` | Distance metric used (Kosman, Euclidean, p-distance, or average squared difference) |

---

## PCAone (Principal Component Analysis)

Principal component analysis (PCA) was performed using PCAone version [VERSION] (Zhang & Meisner 2023) on the unlinked SNP dataset in PLINK format ([N_SNPS] SNPs, [N_SAMPLES] individuals). The [PCAONE_N_PC] leading principal components were computed using SVD method [PCAONE_SVD] (0 = IRAM, 1 = single-pass randomised, 2 = window-based randomised, 3 = full SVD). Eigenvalues and eigenvectors were used to visualise population structure in PC space.

<!-- If PCAone EMU mode was also run: -->
**[Optional EMU variant]** To assess robustness to missing data, PCA was additionally performed using the Expectation–Maximisation (EMU) algorithm as implemented in PCAone, which handles missing genotypes without imputation. The [EMU_N_PC] leading PCs were retained.

<!-- If PCA was run across multiple missing-data thresholds: -->
**[Optional missing-data threshold comparison]** To evaluate sensitivity of PC structure to missing data, the unlinked SNP dataset was filtered to several maximum missing-data thresholds (F_MISSING < [LIST THRESHOLDS from config["parameters"]["PCAone"]["miss"]]) and PCAone was re-run for each filtered subset.

---

## PCoA (Principal Coordinates Analysis)

Principal coordinates analysis (PCoA) was performed on pairwise genetic distance matrices (see [Genetic Distances](05_genetic_distances_networks.md) for distance calculation methods) using the `cmdscale` function in R (R Core Team [VERSION]). The [PCOA_DISTANCE] distance was used. Ordination scores were used to visualise genetic relationships among individuals.

---

## Key software and citations

- **PCAone** – Zhang, C. & Meisner, J. (2023). Fast and accurate out-of-core PCA framework for large biobank data. *Genome Research*, 33, 1599–1608. https://doi.org/10.1101/gr.277525.122
- **R** – R Core Team ([YEAR]). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/
