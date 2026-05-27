# Population Structure and Ancestry Analysis

**Input dataset**: Unlinked SNP dataset (`biallelic_snps_thinned.vcf.gz`; PLINK format `biallelic_snps_thinned.bed/.bim/.fam`; STRUCTURE format `biallelic_snps_thinned.str`)

> If LD pruning or no-thinning strategy was used, replace "unlinked SNP dataset" with "LD-pruned SNP dataset" or "all-biallelic SNP dataset" as appropriate throughout.

---

## Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[N_SNPS]` | SNP count in final analysis VCF |
| `[N_SAMPLES]` | Number of individuals |
| `[K_MIN]`–`[K_MAX]` | `config["parameters"]["k_values"]` (e.g. 1–9) |
| `[REPLICATES]` | `config["parameters"]["structure"]["replicates"]` (default: 10) |
| `[BURNIN]` | `data/mainparams` → `BURNIN` (default: 100,000) |
| `[NUMREPS]` | `data/mainparams` → `NUMREPS` (default: 200,000) |
| `[FS_TOL]` | `config["parameters"]["faststructure"]["tol"]` (default: 10e-6) |
| `[FS_PRIOR]` | `config["parameters"]["faststructure"]["prior"]` (default: `simple`) |
| `[BEST_K_STRUCTURE]` | Optimal K from ΔK / Evanno method (from Evanno plot output) |
| `[BEST_K_ADMIXTURE]` | Optimal K from cross-validation error (from `*.admixture.chooseK_results.txt`) |
| `[BEST_K_FASTSTRUCTURE]` | Optimal K from fastStructure model selection (from `*.chooseK_results.txt`) |
| `[DAPC_N_PCA]` | `config["parameters"]["dapc"]["n_pca"]` (default: 50) |
| `[DAPC_N_DA]` | `config["parameters"]["dapc"]["n_da"]` (default: 10) |
| `[DAPC_CRITERION]` | `config["parameters"]["dapc"]["criterion"]` (default: `diffNgroup`) |
| `[DAPC_BEST_K]` | Optimal cluster number from DAPC cross-validation |
| `[CONSTRUCT_K_MAX]` | Maximum K tested in conStruct |
| `[CONSTRUCT_N_ITER]` | `config["parameters"]["construct"]["n_iterations"]` (default: 10,000) |
| `[CONSTRUCT_N_CHAINS]` | `config["parameters"]["construct"]["n_chains"]` (default: 1) |

---

## STRUCTURE

Bayesian clustering was performed using STRUCTURE version 2.3.4 (Pritchard et al. 2000) to infer the number of genetic clusters (K) and estimate individual ancestry proportions. Analyses were run on the unlinked SNP dataset ([N_SNPS] SNPs, [N_SAMPLES] individuals) using the admixture model with correlated allele frequencies (FREQSCORR = 1, NOADMIX = 0, INFERALPHA = 1), a burn-in of [BURNIN] MCMC iterations, and [NUMREPS] sampling iterations thereafter. K values from [K_MIN] to [K_MAX] were each run [REPLICATES] times with independent random seeds. Run parameters were specified in `data/mainparams` and `data/extraparams`. Replicate runs per K were aligned using the FullSearch algorithm in pophelper v2.3.1 (Francis 2017), and the optimal K was selected using the ΔK method of Evanno et al. (2005), supported by inspection of mean log-likelihood L(K) across runs. The best-supported partition was K = [BEST_K_STRUCTURE].

---

## fastStructure

Ancestry inference was additionally performed using fastStructure version 1.0 (Raj et al. 2014), a variational Bayesian approximation to the STRUCTURE model suited to large SNP datasets. Analyses used the unlinked SNP dataset in PLINK format ([N_SNPS] SNPs, [N_SAMPLES] individuals) with the [FS_PRIOR] prior and a convergence tolerance of [FS_TOL]. K values from [K_MIN] to [K_MAX] were each run once. The optimal K was identified using the `chooseK.py` utility distributed with fastStructure, which selects the K that maximises the marginal likelihood of the data. The best-supported partition was K = [BEST_K_FASTSTRUCTURE].

---

## ADMIXTURE

Maximum-likelihood ancestry estimation was carried out using ADMIXTURE version 1.3.0 (Alexander et al. 2009) on the unlinked SNP dataset in PLINK format ([N_SNPS] SNPs, [N_SAMPLES] individuals). K values from [K_MIN] to [K_MAX] were tested, with cross-validation (CV) error computed using five-fold cross-validation. The optimal K was selected as the value minimising CV error. The best-supported partition was K = [BEST_K_ADMIXTURE].

---

## DAPC (Discriminant Analysis of Principal Components)

Genetic clustering was explored using Discriminant Analysis of Principal Components (DAPC; Jombart et al. 2010) as an unsupervised method independent of any population model. Analyses were performed in R using adegenet version [VERSION] (Jombart 2008; Jombart & Ahmed 2011) on the unlinked SNP dataset ([N_SNPS] SNPs, [N_SAMPLES] individuals). Prior to discriminant analysis, [DAPC_N_PCA] principal component axes were retained. The optimal number of clusters (K) was determined by K-means cross-validation with 30 replicates using the "[DAPC_CRITERION]" criterion across K = [K_MIN]–[K_MAX], identifying K = [DAPC_BEST_K] as the best-supported partition. [DAPC_N_DA] discriminant functions were retained for display.

---

## conStruct (Spatial Population Structure)

<!-- Include this section only if construct analysis was run -->
Spatially explicit population structure was modelled using conStruct version 1.0.4 (Bradburd et al. 2018), which jointly estimates non-spatial and spatial covariance components across K ancestry layers while accounting for isolation by distance. Analyses were conducted on the unlinked SNP dataset ([N_SNPS] SNPs), testing K = 1 to [CONSTRUCT_K_MAX] spatial layers with [CONSTRUCT_N_CHAINS] MCMC chain(s) of [CONSTRUCT_N_ITER] iterations each. Allele frequencies were estimated directly from genotype dosages. Sample coordinates were specified in WGS84 (EPSG:4326).

---

## Key software and citations

- **STRUCTURE** – Pritchard, J.K., Stephens, M. & Donnelly, P. (2000). Inference of population structure using multilocus genotype data. *Genetics*, 155, 945–959. https://doi.org/10.1093/genetics/155.2.945
- **pophelper** – Francis, R.M. (2017). pophelper: an R package and web app to analyse and visualise population structure. *Molecular Ecology Resources*, 17, 27–32. https://doi.org/10.1111/1755-0998.12509
- **Evanno ΔK** – Evanno, G., Regnaut, S. & Goudet, J. (2005). Detecting the number of clusters of individuals using the software STRUCTURE. *Molecular Ecology*, 14, 2611–2620. https://doi.org/10.1111/j.1365-294X.2005.02553.x
- **fastStructure** – Raj, A., Stephens, M. & Pritchard, J.K. (2014). fastSTRUCTURE: variational inference of population structure in large SNP datasets. *Genetics*, 197, 573–589. https://doi.org/10.1534/genetics.114.164350
- **ADMIXTURE** – Alexander, D.H., Novembre, J. & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19, 1655–1664. https://doi.org/10.1101/gr.094052.109
- **adegenet / DAPC** – Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. *Bioinformatics*, 24, 1403–1405; Jombart, T., Devillard, S. & Balloux, F. (2010). Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. *BMC Genetics*, 11, 94. https://doi.org/10.1186/1471-2156-11-94
- **conStruct** – Bradburd, G.S., Coop, G.M. & Ralph, P.L. (2018). Inferring continuous and discrete population genetic structure across space. *Genetics*, 210, 33–52. https://doi.org/10.1534/genetics.118.301333
