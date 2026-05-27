# Genetic Distance Matrices and Networks

**Input dataset**: Unlinked SNP dataset in PLINK format (`biallelic_snps_thinned.bed/.bim/.fam`)

---

## Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[N_SNPS]` | SNP count in final analysis VCF |
| `[N_SAMPLES]` | Number of individuals |
| `[NN_COLOR_BY]` | `config["parameters"]["neighbornet"]["color_by"]` (e.g. Site) |

---

## Genetic distance matrices

Four pairwise individual-level genetic distance metrics were computed from the unlinked SNP dataset ([N_SNPS] SNPs, [N_SAMPLES] individuals) using custom Python scripts operating on PLINK binary files via the bed-reader library:

**Kosman–Leonard distance** (Kosman & Leonard 2005): a standardised pairwise genetic distance for diploid codominant markers. For biallelic dosages (0, 1, 2 copies of the alternate allele), the per-locus contribution is the absolute allele-copy difference divided by two, averaged across loci with data for both individuals.

**Euclidean distance**: pairwise Euclidean distance computed from per-SNP dosage vectors (0, 1, 2), with missing values mean-imputed per locus prior to calculation.

**p-distance** (normalised Hamming distance): for diploid biallelic dosages, per-locus dissimilarity is |g_i − g_j| / 2, averaged across loci where both samples have non-missing data.

**Average squared difference** (bed2diffs-style): pairwise squared genetic differences computed as d_ij = S_ii + S_jj − 2 S_ij, where S = G G^T / n_loci and G is the centred dosage matrix.

---

## NeighborNet split network

A split network was constructed from the p-distance matrix using the NeighborNet algorithm (Bryant & Moulton 2004) as implemented in the phangorn package (Schliep 2011) and visualised using tanggle (Schliep et al. 2023) in R. Tips were coloured by [NN_COLOR_BY].

---

## Key software and citations

- **Kosman–Leonard distance** – Kosman, E. & Leonard, K.J. (2005). Similarity coefficients for molecular markers in studies of genetic relationships between individuals for haploid, diploid, and polyploid species. *Molecular Ecology*, 14, 415–424. https://doi.org/10.1111/j.1365-294X.2004.02416.x
- **NeighborNet** – Bryant, D. & Moulton, V. (2004). Neighbor-Net: An agglomerative method for the construction of phylogenetic networks. *Molecular Biology and Evolution*, 21, 255–265. https://doi.org/10.1093/molbev/msh018
- **phangorn** – Schliep, K.P. (2011). phangorn: phylogenetic analysis in R. *Bioinformatics*, 27, 592–593. https://doi.org/10.1093/bioinformatics/btq706
- **tanggle** – Schliep, K. et al. (2023). tanggle: visualization of split networks. R package. https://github.com/KlausVigo/tanggle
- **bed-reader** – https://github.com/fastlmm/bed-reader
