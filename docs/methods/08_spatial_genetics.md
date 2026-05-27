# Spatial Genetics

---

## MAPI (Mapping Averaged Pairwise Information)

**Input dataset**: Genetic distance matrix (Euclidean distance from unlinked SNP dataset) + sample coordinates

### Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[N_SAMPLES]` | Number of individuals |
| `[N_SNPS]` | SNP count in unlinked SNP dataset |
| `[GRID_HALFWIDTH]` | `config["parameters"]["mapi"]["grid_halfwidth"]` (metres; default: 100,000) |
| `[N_PERM]` | `config["parameters"]["mapi"]["n_permutations"]` (default: 1,000) |
| `[ALPHA]` | `config["parameters"]["mapi"]["alpha"]` (default: 0.05) |
| `[CRS_PROJ]` | `config["parameters"]["mapi"]["crs_projected"]` (default: EPSG:8857 Equal Earth) |

### Methods text

Spatial patterns of genetic differentiation were explored using MAPI (Mapping Averaged Pairwise Information; Rezende et al. 2021). MAPI assigns pairwise genetic distances between individuals to a hexagonal grid of cells covering the study area and computes a cell-level average of all pairwise values for individuals whose midpoint falls within each cell. Input was the Euclidean genetic distance matrix ([N_SAMPLES] × [N_SAMPLES]) derived from the unlinked SNP dataset ([N_SNPS] SNPs). Sample coordinates were projected to [CRS_PROJ] for distance calculations, with a hexagonal grid cell half-width of [GRID_HALFWIDTH] m. Significance of upper and lower tails of the cell-value distribution was assessed by [N_PERM] permutations of individual labels, with cells designated as significant at α = [ALPHA]. Cells with significantly elevated pairwise distances (upper tail) are interpreted as barriers to gene flow, while cells with significantly reduced distances (lower tail) represent corridors of high genetic connectivity.

---

## Key software and citations

- **MAPI** – Rezende, F. et al. (2021). MAPI: an R package for mapping averaged pairwise information. *Methods in Ecology and Evolution*, 12, 1305–1310. https://doi.org/10.1111/2041-210X.13562
