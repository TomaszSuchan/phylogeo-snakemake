# phylogeo-snakemake

Snakemake pipeline for population genetics and phylogeographic data analysis from ipyrad output (`.vcf` and `.loci` files).

## Inputs

- **Main configuration**: `config/config.yaml`
- **Per-project ipyrad output**: `<ipyrad_prefix>.vcf` or `<ipyrad_prefix>.vcf.gz`, and (for pixy) `<ipyrad_prefix>.loci`  
  - `ipyrad_prefix` is defined per project in `config/config.yaml` under `config["projects"][project_name]["ipyrad_prefix"]`
- **Population map (individual → population)**: `data/popmap.txt`
- **Population metadata with coordinates** (optional but required for map-based analyses): `data/popdata.txt`
- **STRUCTURE run parameters** (optional, if `structure` analysis is enabled): `data/mainparams`, `data/extraparams`

## Configuration

All settings are in `config/config.yaml`:
- Define multiple projects under `config["projects"]`, each with its own `config["projects"][project_name]["ipyrad_prefix"]` and `config["projects"][project_name]["analyses"]`
- Global parameters live under `config["parameters"]` and are inherited by all projects unless overridden in `config["projects"][project_name]["parameters"]`
- Enable/disable analyses per project in `config["projects"][project_name]["analyses"]`

Project-specific overrides can be defined only for the settings that differ from the global defaults. For example:

```yaml
parameters:
  vcf_filtering:
    f_missing: 0.5
    maf: 0.0

projects:
  MyProject:
    ipyrad_prefix: "data/myproject"
    analyses:
      pcaone: true
    parameters:
    <<: *global_params
      vcf_filtering:
        f_missing: 0.2
```

In this example, all projects inherit `config["parameters"]["vcf_filtering"]`, but `MyProject` overrides only `config["projects"]["MyProject"]["parameters"]["vcf_filtering"]["f_missing"]`.

## Running

Best is to use conda/mamba:
```
snakemake --use-conda --conda-frontend mamba
```

On HPC you need to specify the executor and resources, eg.:

```
snakemake --executor slurm --cores 10 --default-resource slurm_account=plgapollo2-cpu slurm_partition=plgrid --use-conda --conda-frontend mamba --jobs 10
```

Other useful flags are `--rerun-incomplete`,  `--rerun-triggers mtime`. Add `-n` for a dry-run.

## Workflow

### Preprocessing and Filtering Pipeline

The pipeline applies filtering steps in the following order:

#### 1. Initial VCF Processing
- **Sort VCF**: Sort, compress, and index ipyrad VCF output
- **Index VCF**: Create index file for sorted VCF

#### 2. Sample Subsetting (Optional)
- **User-defined sample subset**: If `config["projects"][project_name]["samples"]` is provided (comma-separated list), only those samples are retained
- **Parameter path**: `config["projects"][project_name]["samples"]`
- **Type**: STRING (comma-separated sample IDs)

#### 3. Relatedness Filtering (Optional)
- **KING-based filtering**: Filter related/clonal individuals using plink2 KING method
- **Parameters** (from `config["parameters"]["relatedness_filtering"]`):
  - `enabled` <BOOLEAN>: Enable/disable relatedness filtering (default: `false`)
  - `king_threshold` <FLOAT>: KING kinship threshold for filtering (default: `0.0884`, equivalent to PI-HAT > 0.2, 2nd-degree relatives)
    - Common values: `0.354` (duplicates/twins), `0.177` (1st-degree), `0.0884` (2nd-degree), `0.0442` (3rd-degree)
  - `mac_threshold` <INTEGER>: Minimum Minor Allele Count (MAC) for KING calculation (default: `1`)
  - `group_by` <LIST[STRING]>: List of columns from popdata to group barplots by (default: `["Region"]`)
- **Output**: Generates filtered sample list and KING relatedness table; creates barplots showing removed individuals grouped by specified columns

#### 4. Missing Data Filtering (After Relatedness)
- **Variant-level missing data filter**: Applied after relatedness filtering (if enabled)
- **Parameters** (from `config["parameters"]["vcf_filtering"]`):
  - `f_missing` <FLOAT>: Maximum missing data threshold for variant filtering (default: `0.5`)
    - Variants missing in > `f_missing` proportion of samples are filtered out
    - `1` keeps all variants, `0` removes variants with any missing data
- **Note**: If `f_missing >= 1`, missing data filtering is skipped

#### 5. Biallelic SNP Filtering
- **Filter**: Retains only biallelic SNPs with Minor Allele Count (MAC) > 1 and, optionally, a minor allele frequency filter
- **Parameters** (from `config["parameters"]["vcf_filtering"]`):
  - `maf` <FLOAT>: Minimum minor allele frequency threshold for variant filtering (default: `0`, i.e. no MAF filter)
    - Implemented as `AF[0] >= maf && AF[0] <= 1-maf` in `bcftools`

#### 6. VCF Thinning (One SNP per RAD Fragment)
- **Thinning / LD-pruning strategy**: Controlled by `config["parameters"]["thinning_strategy"]`:
  - `"thinning"`: Select one SNP per RAD fragment to ensure unlinked markers (default in example config)
  - `"ld_pruning"`: Perform LD-based pruning using `plink` (see LD pruning parameters below)
  - `"none"`: Keep all biallelic SNPs without thinning or pruning
- **Parameters for `"thinning"`** (from `config["parameters"]["vcf_thinning"]`):
  - `min_coverage` <INTEGER>: Minimum number of samples with data at a locus (default: `0`)
  - `method` <STRING>: SNP selection method - `random`, `max_coverage`, or `weighted` (default: `max_coverage`)
  - `ties` <STRING>: How to handle tied coverage values - `first` or `random` (default: `random`)
  - `ns_tag` <STRING>: VCF INFO field tag for number of samples with data (default: `NS`)
  - `id_pattern` <STRING>: Regex pattern to extract RAD fragment ID from variant names (default: `loc(\d+)_`)
- **Parameters for `"ld_pruning"`** (from `config["parameters"]["ld_pruning"]`):
  - `r2` <FLOAT>: \(r^2\) threshold for LD pruning (default: `0.5`)
  - `window_size` <INTEGER>: Window size in kb for LD pruning (default: `50`)
  - `step_size` <INTEGER>: Step size in kb for LD pruning (default: `5`)
- **Key outputs** (per project, under `results/<project>/filtered_data/`):
  - `.../{project}.filtered.vcf.gz`: VCF after missing data filtering
  - `.../{project}.biallelic_snps.vcf.gz`: biallelic SNPs (before thinning / LD pruning)
  - `.../{project}.biallelic_snps_thinned.vcf.gz` or `.../{project}.biallelic_snps_ld_pruned.vcf.gz` (depending on strategy)
  - Corresponding PLINK files: `.../{project}.biallelic_snps_{thinned|ld_pruned|all_snps}.bed/.bim/.fam`

#### 7. Format Conversion
- **PLINK format**: Export thinned VCF to PLINK format (`.bed/.bim/.fam`)
- **Structure format**: Export thinned VCF to Structure format (`.str`)
- **Missing data filtered PLINK files**: Generate PLINK files filtered by missing data thresholds defined in `config["parameters"]["PCAone"]["miss"]` for PCAone robustness analysis
  - Output PLINK prefixes are written under `results/<project>/filtered_data/`

#### Summary of Filtering Order
1. Sort and index input VCF
2. (Optional) Subset to user-specified samples
3. (Optional) Filter related individuals using KING method
4. Filter variants by missing data threshold (`f_missing`)
5. Filter to biallelic SNPs with MAC > 1 (and optional MAF filter)
6. Apply thinning / LD-based pruning strategy (`thinning_strategy`)
7. Export to various formats (PLINK, Structure)
8. (For PCAone) Generate additional VCFs filtered by multiple missing data thresholds

### VCF statistics and missingness

- **Inputs**: raw, filtered and thinned VCFs from the preprocessing pipeline
- **Outputs** (per project, under `results/<project>/stats_vcf/`):
  - `original/`: missingness tables (`*.imiss`, `*.lmiss`), summary stats, and histograms for the original (subset) VCF
  - `filtered/`: the same set of files for the post-filtering, pre-thinning VCF
  - `thinned/`: the same set of files for the thinned / LD-pruned VCF

### Analysis modules and tools descriptions

The sections below summarise all analysis modules currently implemented in the workflow. Each entry states the main software, biological question, core inputs, main outputs, and the most relevant parameters.

For precision, the workflow uses the following internal genotype datasets:

- `raw_sorted.vcf.gz`: the original ipyrad VCF, sorted, compressed, and indexed.
- `subset.vcf.gz`: the raw VCF after optional user-defined sample subsetting.
- `filtered.vcf.gz`: the subset VCF after optional relatedness filtering and variant missing-data filtering.
- `biallelic_snps.vcf.gz`: the filtered VCF restricted to biallelic SNPs with `MAC > 1` and the optional `maf` filter.
- `biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`: the final analysis SNP set after applying `thinning_strategy`. This is what `get_filtered_vcf_output()` returns and is the main input for most downstream SNP-based analyses.
- `biallelic_snps_thinned.bed/.bim/.fam` and `biallelic_snps_thinned.str`: PLINK and STRUCTURE exports created from the final analysis VCF above. The filenames keep the historical `thinned` basename even when the underlying SNP set actually comes from `ld_pruning` or `thinning_strategy: none`.

### Population structure and ancestry

**STRUCTURE** uses Bayesian clustering to estimate ancestry proportions and infer the number of genetic clusters (`K`). It takes the STRUCTURE-format file exported from the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`, written to `biallelic_snps_thinned.str`) together with `data/mainparams` and `data/extraparams`, and produces per-run output files, aligned `Q` matrices, Evanno summaries, barplots, and optional geographic plots. The main parameters are `config["parameters"]["k_values"]` and `config["parameters"]["structure"]["replicates"]`.

**fastStructure** is a faster variational-Bayesian alternative to STRUCTURE for ancestry estimation across multiple `K` values. It uses the PLINK export of the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`, written to `biallelic_snps_thinned.bed/.bim/.fam`) and outputs `.meanQ`, `.meanP`, and `chooseK` summaries, plus downstream barplots and maps when metadata are available. The main parameters are `config["parameters"]["k_values"]`, `config["parameters"]["faststructure"]["tol"]`, and `config["parameters"]["faststructure"]["prior"]`.

**ADMIXTURE** estimates ancestry proportions by maximum likelihood and reports cross-validation error for model choice across `K`. It uses the PLINK export of the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`, written to `biallelic_snps_thinned.bed/.bim/.fam`) and produces `.Q` and `.P` files together with downstream barplots, maps, and `K`-selection summaries. The main parameter is `config["parameters"]["k_values"]`.

**evalAdmix** evaluates the fit of ancestry models by calculating correlations of residuals after ADMIXTURE, fastStructure, or STRUCTURE inference. It uses the same PLINK export of the final analysis VCF together with the corresponding `Q` and `P` matrices and produces residual-correlation files (`.corres`) and diagnostic plots under the relevant ancestry-method directory. This is especially useful for checking whether the inferred admixture model adequately captures covariance among individuals.

**mapmixture** is the shared geographic visualisation layer for ancestry analyses. It uses ancestry coefficients from STRUCTURE, fastStructure, ADMIXTURE, and DAPC together with `indpopdata.txt` to generate population pie-chart maps and barplots. The main parameters are in `config["parameters"]["mapmixture"]`, with shared map settings in `config["parameters"]["map_background"]`.

**conStruct** models spatial population structure by allowing covariance to decay continuously with geography rather than assuming strictly discrete clusters. It uses the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`) and `indpopdata.txt`, and produces per-`K` `.rds` result files, layer-proportion tables, and map/barplot visualisations. The main parameters are `config["parameters"]["k_values"]` and the block `config["parameters"]["construct"]` (`n_chains`, `n_iterations`, `make.freqs`, `geoDist`, `coords`, `save.files`, `crs`).

**DAPC** (Discriminant Analysis of Principal Components) provides multivariate clustering and assignment based on the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`). It produces criterion/BIC summaries, cluster assignments, membership probabilities, scatterplots, and optional geographic maps. The main parameters are `config["parameters"]["k_values"]` and the block `config["parameters"]["dapc"]` (`n_pca`, `n_da`, `criterion`).

### Ordination and genetic distances

**PCAone** performs principal component analysis on the PLINK export of the final analysis VCF and is intended as the main fast PCA implementation for large SNP matrices. It produces eigenvector/eigenvalue files and a consistent set of single-panel and faceted plots. The main parameters are `config["parameters"]["PCAone"]["SVD_method"]` and `config["parameters"]["PCAone"]["PCnum"]`, while plotting is controlled by `config["parameters"]["pca_plot"]` (`pc_max`, `color_by`, `include_missing`, `pca_colors`).

**PCAone with EMU** uses the EMU mode within PCAone to model missing genotypes during PCA. It uses the same inputs and output structure as standard PCAone, but is preferable when missingness is substantial and should be handled during ordination rather than by prior filtering. It uses the same `config["parameters"]["PCAone"]` and `config["parameters"]["pca_plot"]` settings as PCAone.

**PCAone with missing-data thresholds** repeats PCAone after creating separate SNP datasets filtered at each value in `config["parameters"]["PCAone"]["miss"]`. This allows explicit sensitivity analysis of ordination results to marker missingness. Each threshold produces separate PLINK files, eigenvectors, eigenvalues, and plot sets.

**EMU-PCA** runs the standalone `emu` program rather than PCAone's internal EMU mode. It is intended for low-coverage or uncertainty-prone genotype data where expectation-maximisation-based PCA is preferred. It uses the PLINK export of the final analysis VCF and produces EMU eigenvectors/eigenvalues plus the same plotting family. The main parameters are `config["parameters"]["EMU"]` (`eig`, `eig_out`) and `config["parameters"]["pca_plot"]`.

**VCF2PCACluster** performs PCA with its own internal SNP filtering and optional clustering/kinship correction, using the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`) as input. It therefore starts from the main downstream SNP set and then applies its own internal `MAF`, missingness, heterozygosity, and HWE filters. It outputs eigenvectors and eigenvalues for every tested `MAF × Miss` combination. The main parameters are in `config["parameters"]["vcf2pcacluster"]`, especially `cluster_method`, `PCnum`, `KinshipMethod`, and the nested block `config["parameters"]["vcf2pcacluster"]["SNP_filtering"]` (`MAF`, `Miss`, `Het`, `HWE`, `Fchr`).

**Euclidean genetic distance** computes pairwise distances between individuals from the PLINK export of the final analysis VCF by treating biallelic genotypes as alternate-allele dosages (`0`, `1`, `2`) and mean-imputing missing values per SNP. It is implemented in Python with `bed-reader` and `scipy` and produces `results/<project>/gen_dist/<project>.euclidean_distance.tsv`. This matrix is the downstream input for PCoA and MAPI.

**PCoA** performs principal coordinates analysis on the Euclidean distance matrix. It outputs eigenvector/eigenvalue files in a PCA-like format and reuses the same plotting system as PCAone for colored, labeled, missingness-based, and faceted ordinations. Plot behaviour is controlled by `config["parameters"]["pca_plot"]`.

### Genetic diversity, relatedness, and demographic summaries

**Pixy** calculates unbiased estimates of nucleotide diversity (`pi`), genetic differentiation (`FST`), and sequence divergence (`dXY`) using both variant and invariant sites extracted from the original ipyrad `.loci` file rather than from the downstream SNP-only VCFs. It uses the project keep-list after sample filtering plus grouping information from `indpopdata.txt`, and produces raw statistic tables, bootstrap summaries, heatmaps, barplots, and optional diversity maps. The main parameters are in `config["parameters"]["pixy"]` (`stats`, `window_size`, `bootstrap_replicates`, `group_by`, `point_size`, `map_outline`).

**Genome scan** performs sliding-window analyses of `FST`, `pi`, and `dXY` between two user-defined groups. It uses an invariant-site VCF reconstructed from the original ipyrad data, together with `indpopdata.txt`, to identify samples whose `config["parameters"]["genome_scan"]["population_column"]` value matches `pop1` or `pop2`. The workflow then subsets the invariant-site VCF to retain only individuals from those two groups, creates a two-population popmap from the same subset, and runs `pixy` on that reduced dataset; it therefore does not use the final thinned SNP set and does not include samples outside the focal comparison. It produces statistic tables and genome-wide plots under `results/<project>/genome_scan/`. The main parameters are in `config["parameters"]["genome_scan"]` (`population_column`, `pop1`, `pop2`, `window_size_fst`, `window_size_pi`, `window_size_dxy`).

**Relatedness** estimates pairwise kinship using multiple estimators, including the Yang et al. genomic relatedness statistic, Manichaikul/KING-style kinship, PLINK IBD summaries, and PLINK2 KING outputs. The main relatedness analysis uses the final analysis VCF and its derived PLINK files, whereas the earlier relatedness-filtering step in preprocessing uses `subset.vcf.gz` before variant missingness filtering to decide which individuals to retain. Outputs are `.relatedness`, `.relatedness2`, `.genome`, and `.king` tables, plus metadata-coloured network plots. Visualisation is controlled by `config["parameters"]["relatedness_plot"]` (`color_by`, `relatedness_colors`).

**ROH** uses `bcftools roh` to identify runs of homozygosity and summarises them across individuals and groups. It uses the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`) together with `indpopdata.txt`, and produces raw ROH calls, per-sample summaries, comparison tables, and grouped plots. The main parameter is `config["parameters"]["roh"]["group_by"]`.

**AMOVA** partitions molecular variance across hierarchical levels defined in the metadata. It uses the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`) and `indpopdata.txt`, and produces an AMOVA results table and variance-component plot. The main parameters are `config["parameters"]["amova"]["strata"]` and `config["parameters"]["amova"]["nperm"]`.

### Spatial landscape genetics

**MAPI** (Mapping Averaged Pairwise Information) identifies spatial regions where observed genetic distances are unusually high or low relative to geography. It uses the Euclidean distance matrix and `indpopdata.txt`, and produces result and significance-tail GeoPackages plus a final PDF map. The main parameters are in `config["parameters"]["mapi"]` (`n_permutations`, `grid_halfwidth`, `crs_projected`, `crs_geographic`, `alpha`, `fill_var`).

### Phylogenetic and historical inference

**IQ-TREE** performs maximum-likelihood phylogenetic inference from the ipyrad `.phy` alignment. It outputs tree files, run logs, support summaries, and plotted trees. The main parameters are in `config["parameters"]["iqtree"]` (`model`, `bootstraps`, optional `outgroup`, `robust-phy`, `support_threshold`).

**trimAl + IQ-TREE** adds an alignment-trimming step before phylogenetic inference. It uses `trimAl` to remove sites with excessive gaps and then runs IQ-TREE on the trimmed alignment, producing trimmed alignments together with the same tree and plotting outputs as the untrimmed workflow. The main trimming parameter is `config["parameters"]["trimal"]["gt"]`; the downstream tree step uses the same `config["parameters"]["iqtree"]` settings.

**fineRADstructure** uses `RADpainter` and `fineRADstructure` to infer recent coancestry and hierarchical relationships among individuals from RAD-seq loci. It takes `biallelic_snps.vcf.gz`, that is, the biallelic SNP set before thinning or LD pruning, together with `indpopdata.txt`, and produces the fineRADstructure input matrix, coancestry chunks, MCMC and tree XML outputs, and individual/population-level coancestry heatmaps. The main parameters are in `config["parameters"]["fineradstructure"]`, especially the nested `cluster`, `tree`, and `plot` blocks.

### Population metadata

**Population assignment** is generated automatically from either a user-provided `popmap`, a richer `popdata` table with geographic coordinates and grouping variables, or a `popseparator` that parses population codes from sample names. The main inputs are `config["parameters"]["popmap"]`, `config["parameters"]["popdata"]`, and `config["parameters"]["popseparator"]`, and the main outputs are `results/<project>/indpopdata.txt`, `results/<project>/indpopdata_all.txt`, and `results/<project>/stats_samples/<project>.population_summary.txt`.

**Population maps and summaries** use `indpopdata.txt` and the shared map settings to create descriptive sampling maps. If `popdata` is available, the pipeline can generate `results/<project>/stats_samples/<project>.population_map.pdf` together with tabular summaries of sample counts per population. The main map controls are `config["parameters"]["population_map"]`, `config["parameters"]["map_background"]`, and `config["parameters"]["map_boundary"]`.

### Computational resources

Resources per analysis type defined under `config["parameters"]["resources"]`:
- `threads` <INTEGER>: Number of CPU cores
- `mem_mb` <INTEGER>: Memory allocation in megabytes
- `runtime` <INTEGER>: Maximum runtime in minutes