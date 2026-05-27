# Phylogenetic Analysis

---

## IQ-TREE (Maximum-Likelihood Phylogeny)

**Input dataset**: Concatenated phylip alignment from ipyrad (`*.phy`) — all variable sites, no LD filtering

### Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[N_SAMPLES]` | Number of individuals in the alignment |
| `[N_SITES]` | Alignment length (number of sites) |
| `[MODEL]` | `config["parameters"]["iqtree"]["model"]` (default: `MFP` = ModelFinder Plus) |
| `[BOOTSTRAP_N]` | `config["parameters"]["iqtree"]["bootstraps"]` (default: 1,000) |
| `[BEST_MODEL]` | Substitution model selected by ModelFinder (from `*.iqtree` log) |
| `[OUTGROUP]` | `config["parameters"]["iqtree"]["outgroup"]` (if specified; otherwise midpoint rooting) |
| `[SUPPORT_THRESHOLD]` | `config["parameters"]["iqtree"]["support_threshold"]` (default: 70) |

### Methods text

A maximum-likelihood phylogenetic tree was inferred using IQ-TREE version [VERSION] (Minh et al. 2020) from the concatenated RAD-seq alignment produced by ipyrad ([N_SAMPLES] individuals, [N_SITES] sites). The best-fit substitution model was identified by ModelFinder (Kalyaanamoorthy et al. 2017) from the full model space (`-m MFP`) and was [BEST_MODEL]. Branch support was assessed with [BOOTSTRAP_N] ultrafast bootstrap replicates (Hoang et al. 2018). The resulting tree was [rooted on [OUTGROUP] / midpoint-rooted] for display; nodes with bootstrap support below [SUPPORT_THRESHOLD]% are not shown.

---

## trimAl + IQ-TREE

<!-- Include this section only if iqtree_trimal analysis was run -->

### Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[TRIMAL_GT]` | `config["parameters"]["trimal"]["gt"]` (default: 0.95) |
| `[N_SITES_TRIMMED]` | Alignment length after trimming |

### Methods text

Prior to phylogenetic inference, the concatenated alignment was trimmed using trimAl version [VERSION] (Capella-Gutiérrez et al. 2009) with a gap threshold of [TRIMAL_GT] (sites with more than [100*(1-TRIMAL_GT)]% gaps were removed), retaining [N_SITES_TRIMMED] sites. The trimmed alignment was then analysed with IQ-TREE as described above.

---

## fineRADstructure (Coancestry-Based Clustering)

**Input dataset**: Biallelic SNP VCF before LD thinning (`biallelic_snps.vcf.gz`) — uses haplotype information within RAD loci

### Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[N_SAMPLES]` | Number of individuals |
| `[N_CHUNKS]` | Number of coancestry chunks (from RADpainter output) |
| `[MCMC_ITER_CLUSTER]` | `config["parameters"]["fineradstructure"]["cluster"]["mcmc_iterations"]` (default: 1,000,000) |
| `[BURNIN_CLUSTER]` | `config["parameters"]["fineradstructure"]["cluster"]["burnin"]` (default: 1,000,000) |
| `[THIN_CLUSTER]` | `config["parameters"]["fineradstructure"]["cluster"]["thinning"]` (default: 1,000) |
| `[MCMC_ITER_TREE]` | `config["parameters"]["fineradstructure"]["tree"]["mcmc_iterations"]` (default: 10,000) |

### Methods text

Co-ancestry-based population clustering was performed using fineRADstructure (Malinsky et al. 2018). First, RADpainter was used to compute a pairwise haplotype co-ancestry matrix from the biallelic SNP dataset ([N_SAMPLES] individuals, [N_CHUNKS] informative co-ancestry chunks), exploiting haplotype phase information within individual RAD loci. The co-ancestry matrix was then analysed using fineSTRUCTURE (Lawson et al. 2012) via a Markov chain Monte Carlo (MCMC) clustering algorithm with [MCMC_ITER_CLUSTER] sampling iterations ([BURNIN_CLUSTER] burn-in, sampled every [THIN_CLUSTER] iterations). A maximum-clade-credibility population tree was inferred from the MCMC output using [MCMC_ITER_TREE] tree-building iterations. Results were visualised as a co-ancestry heatmap with hierarchical clustering.

---

## Key software and citations

- **IQ-TREE 2** – Minh, B.Q. et al. (2020). IQ-TREE 2: new models and methods for phylogenetic inference. *Molecular Biology and Evolution*, 37, 1530–1534. https://doi.org/10.1093/molbev/msaa015
- **ModelFinder** – Kalyaanamoorthy, S. et al. (2017). ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods*, 14, 587–589. https://doi.org/10.1038/nmeth.4285
- **Ultrafast bootstrap** – Hoang, D.T. et al. (2018). UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution*, 35, 518–522. https://doi.org/10.1093/molbev/msx281
- **trimAl** – Capella-Gutiérrez, S., Silla-Martínez, J.M. & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic studies. *Bioinformatics*, 25, 1972–1973. https://doi.org/10.1093/bioinformatics/btp348
- **fineRADstructure / RADpainter** – Malinsky, M. et al. (2018). RADpainter and fineRADstructure: population inference from RADseq data. *Molecular Biology and Evolution*, 35, 1284–1290. https://doi.org/10.1093/molbev/msy023
- **fineSTRUCTURE** – Lawson, D.J. et al. (2012). Inference of population structure using dense haplotype data. *PLOS Genetics*, 8, e1002453. https://doi.org/10.1371/journal.pgen.1002453
