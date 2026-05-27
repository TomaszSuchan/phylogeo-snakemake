# Genetic Diversity and Differentiation

---

## pixy (π, FST, dXY)

**Input dataset**: All-sites dataset (`allsites.vcf.gz`) — *invariant sites are required for unbiased estimates*

### Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[N_SAMPLES]` | Number of individuals |
| `[N_POPS]` | Number of populations / grouping levels used |
| `[WINDOW_SIZE]` | `config["parameters"]["pixy"]["window_size"]` (default: 10,000 bp) |
| `[BOOTSTRAP_N]` | `config["parameters"]["pixy"]["bootstrap_replicates"]` (default: 1,000) |
| `[GROUP_BY]` | `config["parameters"]["pixy"]["group_by"]` (e.g. Region, Site) |
| `[STATS]` | `config["parameters"]["pixy"]["stats"]` (e.g. pi, fst, dxy) |

### Methods text

Nucleotide diversity (π), pairwise FST, and absolute pairwise divergence (dXY) were estimated using pixy version [VERSION] (Korunes & Samuk 2021) to avoid the ascertainment biases that arise when statistics are computed exclusively from variant sites. Pixy operates on an all-sites VCF that includes both variant and invariant positions, reconstructed from the ipyrad `.loci` file (see [Filtering](01_filtering.md)). Populations were defined by [GROUP_BY] assignments. Statistics were computed in non-overlapping windows of [WINDOW_SIZE] bp. Population-level summaries and confidence intervals were obtained by bootstrapping across windows ([BOOTSTRAP_N] replicates). [FST values are reported as weighted averages across windows (Weir & Cockerham 1984).]

<!-- Include dXY sentence only if dxy was computed -->
Absolute divergence (dXY) between population pairs was calculated in the same window framework and summarised as the genome-wide mean ± 95% bootstrap CI.

---

## AMOVA (Analysis of Molecular Variance)

**Input dataset**: Unlinked SNP dataset (`biallelic_snps_thinned.vcf.gz`)

### Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[STRATA]` | `config["parameters"]["amova"]["strata"]` (e.g. Region, Site) |
| `[NPERM]` | `config["parameters"]["amova"]["nperm"]` (default: 999) |

### Methods text

Analysis of molecular variance (AMOVA; Excoffier et al. 1992) was performed in R using the poppr package (Kamvar et al. 2014) to partition genetic variance among and within hierarchical groupings ([STRATA]). The unlinked SNP dataset was used. Statistical significance of variance components was assessed by [NPERM] random permutations of individuals among strata.

---

## Key software and citations

- **pixy** – Korunes, K.L. & Samuk, K. (2021). pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. *Molecular Ecology Resources*, 21, 1359–1368. https://doi.org/10.1111/1755-0998.13326
- **AMOVA** – Excoffier, L., Smouse, P.E. & Quattro, J.M. (1992). Analysis of molecular variance inferred from metric distances among DNA haplotypes. *Genetics*, 131, 479–491.
- **poppr** – Kamvar, Z.N., Tabima, J.F. & Grünwald, N.J. (2014). Poppr: an R package for genetic analysis of populations with clonal, partially clonal, and/or sexual reproduction. *PeerJ*, 2, e281. https://doi.org/10.7717/peerj.281
