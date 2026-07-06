# phylogeo-snakemake

Snakemake pipeline for population genetics and phylogeographic data analysis from ipyrad output (`.vcf` and `.loci` files).

### How software is installed

Most tools run in per-rule conda environments defined in `workflow/envs/*.yaml`. With `--use-conda`, Snakemake creates those environments in `.snakemake/conda/` (hashed by the YAML contents) before jobs are submitted. Environments use **conda-forge** and **bioconda** only.

A few R packages are not on conda and are installed once by dedicated rules from CRAN or GitHub (for example `mapmixture`, `elevatr`, `alstructure`, `tess3r`, `conStruct`). Those rules run inside the conda env from the matching YAML, then call `install.packages()` or `remotes::install_github()` as needed.

Always use **mamba** as the conda frontend (`--conda-frontend mamba`): it solves environments faster and with less memory than plain `conda`.

**HPC note:** Snakemake builds missing conda environments on the machine that runs the Snakemake driver (often the login node). Some solves (for example `pixy`) can hit login-node memory limits and fail at `Solving environment...` with no clear error. If that happens, create environments on a compute node first, then launch the workflow from the login node as usual:

```bash
# 1) On a compute node with enough RAM and walltime
srun --partition=plgrid --account=your-cpu-account --mem=32G --cpus-per-task=4 --time=120:00 bash
cd phylogeo-snakemake
source /path/to/miniforge3/bin/activate snakemake   # or your snakemake env

snakemake --conda-create-envs-only \
  --use-conda --conda-frontend mamba \
  --configfile config/your_project.yaml

# 2) Back on the login node — envs already exist
exit

# 3) run the snakemake pipeline as usula from the login node

snakemake --executor slurm \
  --default-resources slurm_account=your-cpu-account slurm_partition=plgrid runtime=300 \
  --cores 48 --jobs 300 \
  --use-conda --conda-frontend mamba \
  --configfile config/your_project.yaml
```

On scratch filesystems with inode quotas, run `conda clean -a -y` periodically; conda package caches and many Snakemake envs each hold thousands of files.

## Cloning

This repository uses a git submodule. Clone with:

```bash
git clone --recurse-submodules https://github.com/TomaszSuchan/phylogeo-snakemake.git
```

If you already have a clone without the submodule:

```bash
git submodule update --init --recursive
```

## Inputs

- **Main configuration**: pass explicitly with `--configfile` (e.g. `config/config.yaml`)
- **Per-project ipyrad output**: `<ipyrad_prefix>.vcf` or `<ipyrad_prefix>.vcf.gz`, and (for pixy and genome_scan) `<ipyrad_prefix>.loci`  
  - `ipyrad_prefix` is defined per project in `config/config.yaml` under `projects.<project_name>.ipyrad_prefix`
- **Population map (individual → population)**: `data/popmap.txt`
- **Population metadata with coordinates** (optional but required for map-based analyses): `data/popdata.txt`
- **STRUCTURE run parameters** (optional, if `structure` analysis is enabled): `data/mainparams`, `data/extraparams`

## Configuration

All settings are in `config/config.yaml`. Global defaults live under `parameters` and are inherited by every project. Project-specific paths, analysis flags, and overrides live under `projects.<project_name>`.


**Config path notation:** In this README, settings are cited with dot paths that mirror `config/config.yaml`. Global defaults use `parameters.<section>.<key>` (for example `parameters.tess3.ploidy`). Project-specific overrides use `projects.<project_name>.<key>` (for example `projects.MyProject.samples`). Analysis enable flags live under `projects.<project_name>.analyses.<method>`. In Snakefile code the same values are read as `config["parameters"][...]` or `config["projects"][...]`.

Start by copying the commented `my_study` template in `config/config.yaml`, then edit only the values that differ from the defaults. A minimal project looks like this:

```yaml
projects:
  MyProject:
    ipyrad_prefix: "data/myproject"
    analyses:
      pcaone: true
      pixy: false
      structure: false
    parameters:
      <<: *global_params
      popmap: "data/myproject.popmap.txt"
      popdata: "data/myproject.popdata.txt"
      vcf_filtering:
        f_missing: 0.2
```

In this example, `MyProject` inherits every global parameter, then overrides only its input metadata paths and missing-data threshold.

### Choosing Project Inputs

- `projects.<project_name>.ipyrad_prefix`: path prefix before the extension. If this is `data/myproject`, the workflow looks for `data/myproject.vcf` or `data/myproject.vcf.gz`; pixy and genome-scan analyses also need `data/myproject.loci`.
- `projects.<project_name>.samples`: optional comma-separated sample IDs used to subset the whole project before filtering. Omit it or leave it empty to use all samples in the VCF.
- `parameters.popmap`: two-column individual-to-population file. Use this when you only need population assignment.
- `parameters.popseparator`: separator used to parse population names from sample IDs when no `popmap` is available.
- `parameters.popdata`: richer metadata table used for maps, grouping, colors, labels, and stratified summaries. Required columns are `Site`, `Lat`, and `Lon`; additional columns can be used by parameters such as `parameters.pca_plot.color_by`, `parameters.pixy.group_by`, `parameters.pixy.arrange_by`, `parameters.mapmixture.population_labels`, and `parameters.amova.strata`.

Column names are matched exactly, including spelling and case. If `popdata` contains `Region`, use `Region`, not `region`.

### Choosing Analyses

Per-project analysis flags include `structure`, `faststructure`, `admixture`, `evaladmix`, `tess3`, `snmf`, `alstructure`, `construct`, `treemix`, `dapc`, `spca`, `pcaone`, `pcaone_emu`, `pcaone_miss`, `emupca`, `vcf2pcacluster`, `gen_dist`, `pcoa`, `neighbornet`, `mapi`, `pixy`, `genome_scan`, `relatedness`, `roh`, `amova`, `iqtree`, `iqtree_trimal`, `iqtree_robust`, `njtree`, and `fineradstructure`. Set a flag to `true` to request that module's final outputs.

Some flags need companion settings:

- `evaladmix` needs at least one ancestry method also enabled: `admixture`, `faststructure`, or `structure`.
- `tess3`, `construct`, `dapc` maps, `spca`, `mapi`, population maps, ancestry maps, pixy maps, and other geographic outputs need `popdata` with usable coordinates.
- `gen_dist` is required for `pcoa`, `mapi`, and NeighborNet distance inputs; the workflow requests needed distance files automatically when those modules are enabled.
- `pixy` and `genome_scan` need the original ipyrad `.loci` file because they reconstruct invariant sites.
- `njtree` needs `{ipyrad_prefix}.phy` in the ipyrad output directory (the concatenated PHYLIP alignment exported by ipyrad).
- `genome_scan` needs `parameters.genome_scan.population_column`, `parameters.genome_scan.pop1`, and `parameters.genome_scan.pop2`.
- `amova` needs a non-empty `parameters.amova.strata` list.
- `treemix` groups samples by `parameters.treemix.population_column`, which must be a column in `indpopdata.txt` such as `Site` or a metadata column from `popdata`.

### Choosing Resources

Resources are configured under `parameters.resources`. Each block uses:

- `parameters.resources.<module>.threads` (and similarly `mem_mb`, `runtime`). Example fields:

- `parameters.resources.<module>.threads`: CPU cores requested by the rule.
- `parameters.resources.<module>.mem_mb`: memory in megabytes.
- `parameters.resources.<module>.runtime`: walltime in minutes, especially useful for cluster executors.

Override resources per project when your dataset is larger or your cluster has different limits. Expensive modules commonly needing adjustment include `structure`, `construct`, `pixy`, `pcaone`, `iqtree`, `fineradstructure`, and `genome_scan`.

## Running

Use conda/mamba for all rules (`--use-conda --conda-frontend mamba`). See **How software is installed** above if environment creation fails on the login node.

Local or interactive run:

```bash
snakemake --configfile config/config.yaml --use-conda --conda-frontend mamba
```

On HPC, specify the executor and resources:

```bash
snakemake --configfile config/config.yaml --executor slurm --cores 10 \
  --default-resources slurm_account=plgapollo2-cpu slurm_partition=plgrid \
  --use-conda --conda-frontend mamba --jobs 10
```

Other useful flags: `--rerun-incomplete`, `--rerun-triggers mtime`. Add `-n` for a dry-run.

Optional: set strict channel priority once in your conda installation:

```bash
conda config --set channel_priority strict
```

## Workflow

### Preprocessing and Filtering Pipeline

The pipeline applies filtering steps in the following order:

#### 1. Initial VCF Processing
- **Sort VCF**: Sort, compress, and index ipyrad VCF output
- **Index VCF**: Create index file for sorted VCF

#### 2. Sample Subsetting (Optional)
- **User-defined sample subset**: If `projects.<project_name>.samples` is provided (comma-separated list), only those samples are retained. Use this to run the whole workflow on selected individuals; because it happens before all filters, downstream summaries only see the retained samples.
- **Parameter path**: `projects.<project_name>.samples`
- **Type**: STRING (comma-separated sample IDs)

#### 3. Relatedness Filtering (Optional)
- **KING-based filtering**: Filter related/clonal individuals using plink2 KING method
- **Parameters** (from `parameters.relatedness_filtering`):
  - `parameters.relatedness_filtering.enabled` <BOOLEAN>: Enable/disable relatedness filtering (default: `false`). Enable this when close kin could bias PCA, ancestry clustering, allele frequencies, or diversity summaries.
  - `parameters.relatedness_filtering.king_threshold` <FLOAT>: KING kinship threshold for filtering (default: `0.0884`, equivalent to PI-HAT > 0.2, 2nd-degree relatives)
    - Common values: `0.354` (duplicates/twins), `0.177` (1st-degree), `0.0884` (2nd-degree), `0.0442` (3rd-degree). Lower values are stricter and remove more distant relatives.
  - `parameters.relatedness_filtering.mac_threshold` <INTEGER>: Minimum Minor Allele Count (MAC) for KING calculation (default: `1`). Raising this excludes very rare variants from the kinship calculation, which can help with sparse RAD datasets.
  - `parameters.relatedness_filtering.group_by` <LIST[STRING]>: List of columns from popdata to group barplots by (default: `["Region"]`). These columns must exist in `indpopdata`.
- **Output**: Generates filtered sample list and KING relatedness table; creates barplots showing removed individuals grouped by specified columns

#### 4. Missing Data Filtering (After Relatedness)
- **Variant-level missing data filter**: Applied after relatedness filtering (if enabled)
- **Parameters** (from `parameters.vcf_filtering`):
  - `parameters.vcf_filtering.f_missing` <FLOAT>: Maximum missing data threshold for variant filtering (default: `1.0`)
    - Variants missing in > `f_missing` proportion of samples are filtered out
    - `1` keeps all variants, `0.5` keeps variants genotyped in at least half of samples, and `0` removes variants with any missing data. Lower values are stricter and reduce missingness but also reduce SNP count. This is a site-level filter that drops poorly-genotyped loci, not rare alleles. 
- **Note**: If `f_missing >= 1`, missing data filtering is skipped

#### 5. Biallelic SNP Filtering
- **Filter**: Retains only biallelic SNPs with Minor Allele Count (MAC) > 1 and, optionally, a minor allele frequency filter
- **Parameters** (from `parameters.vcf_filtering`):
  - `parameters.vcf_filtering.maf` <FLOAT>: Minimum minor allele frequency threshold for variant filtering (default: `0`, i.e. no MAF filter)
    - Implemented as `AF[0] >= maf && AF[0] <= 1-maf` in `bcftools`
    - `0` applies no MAF filter, `0.01` is a light filter, and `0.05` focuses on common variants. Higher values can stabilise ordination and clustering but remove rare/local variants. Note that `0` does **not** keep the very rarest variants: the hardcoded `MAC > 1` filter in this step always drops singletons (alleles seen in a single copy) regardless of `maf`.

#### 6. VCF Thinning (One SNP per RAD Fragment)
- **Thinning / LD-pruning strategy**: Controlled by `parameters.thinning_strategy`:
  - `"thinning"`: Select one SNP per RAD fragment to ensure unlinked markers (default)
  - `"ld_pruning"`: Perform LD-based pruning using `plink` (see LD pruning parameters below)
  - `"none"`: Keep all biallelic SNPs without thinning or pruning
  - Most population structure workflows should use `"thinning"` or `"ld_pruning"` to reduce correlations among SNPs.
- **Parameters for `"thinning"`** (from `parameters.vcf_thinning`):
  - `parameters.vcf_thinning.min_coverage` <INTEGER>: Minimum number of samples with data at a locus (default: `0`). Increase it when selected SNPs should be broadly represented across samples.
  - `parameters.vcf_thinning.method` <STRING>: SNP selection method - `random`, `max_coverage`, or `weighted` (default: `max_coverage`). `max_coverage` picks the best-covered SNP per RAD fragment; `weighted` favours better-covered SNPs without always taking the maximum.
  - `parameters.vcf_thinning.ties` <STRING>: How to handle tied coverage values - `first` or `random` (default: `random`). `random` avoids systematic first-SNP bias; `first` is reproducible if input order is stable.
  - `parameters.vcf_thinning.ns_tag` <STRING>: VCF INFO field tag for number of samples with data (default: `NS`)
  - `parameters.vcf_thinning.id_pattern` <STRING>: Regex pattern to extract RAD fragment ID from variant names (default: `loc(\d+)_`). Change it if your VCF IDs do not use ipyrad-style locus IDs.
- **Parameters for `"ld_pruning"`** (from `parameters.ld_pruning`):
  - `parameters.ld_pruning.r2` <FLOAT>: \(r^2\) threshold for LD pruning (default: `0.5`). Lower values are stricter and keep fewer, more independent SNPs.
  - `parameters.ld_pruning.window_size` <INTEGER>: Window size in kb for LD pruning (default: `50`). Larger windows detect broader LD blocks but may remove more SNPs.
  - `parameters.ld_pruning.step_size` <INTEGER>: Step size in kb for LD pruning (default: `5`). This controls how far the window advances between LD checks.
- **Key outputs** (per project, under `results/<project>/filtered_data/`):
  - `.../{project}.filtered.vcf.gz`: VCF after missing data filtering
  - `.../{project}.biallelic_snps.vcf.gz`: biallelic SNPs (before thinning / LD pruning)
  - `.../{project}.biallelic_snps_thinned.vcf.gz` or `.../{project}.biallelic_snps_ld_pruned.vcf.gz` (depending on strategy)
  - Corresponding PLINK files: `.../{project}.biallelic_snps_{thinned|ld_pruned|all_snps}.bed/.bim/.fam`

#### 7. Format Conversion
- **PLINK format**: Export thinned VCF to PLINK format (`.bed/.bim/.fam`)
- **Structure format**: Export thinned VCF to Structure format (`.str`)
- **Missing data filtered PLINK files**: Generate PLINK files filtered by missing data thresholds defined in `parameters.PCAone.miss` for PCAone robustness analysis
  - Output PLINK prefixes are written under `results/<project>/filtered_data/`

#### Summary of Filtering Order
1. Sort and index input VCF
2. (Optional) Subset to user-specified samples
3. (Optional) Filter related individuals using KING method
4. Filter variants by missing data threshold (`parameters.vcf_filtering.f_missing`)
5. Filter to biallelic SNPs with MAC > 1 (and optional MAF filter)
6. Apply thinning / LD-based pruning strategy (`parameters.thinning_strategy`)
7. Export to various formats (PLINK, Structure)
8. (For PCAone) Generate additional VCFs filtered by multiple missing data thresholds

### VCF statistics, missingness, and sequencing depth

- **Inputs**: raw, filtered and thinned VCFs from the preprocessing pipeline
- **Outputs** (per project, under `results/<project>/stats_vcf/`):
  - `original/`: missingness tables (`*.imiss`, `*.lmiss`), sequencing-depth tables, summary stats, and histograms for the original (subset) VCF
  - `filtered/`: the same set of files for the post-filtering, pre-thinning VCF
  - `thinned/`: the same set of files for the thinned / LD-pruned VCF
- **Sequencing depth outputs**:
  - `*.idepth`: per-individual mean depth from `vcftools --depth`, useful for identifying low-coverage samples
  - `*.ldepth.mean`: per-site mean depth from `vcftools --site-mean-depth`, useful for identifying poorly covered or unusually high-depth loci
  - `*.depth_summary.txt`: combined report with overall mean depth across called genotypes, per-individual and per-site mean-depth summaries, and counts below common depth thresholds (`<3`, `<5`, `<10`) or above high depth (`>100`)
- **Reporting guidance**: use the overall mean depth across called genotypes as the headline sequencing-depth value in methods, and inspect per-individual and per-site summaries separately because low-depth samples and low-depth loci have different biological and QC implications.

### Analysis modules and tools descriptions

The sections below summarise all analysis modules currently implemented in the workflow. Each entry focuses on what the tool does, what data it uses, and which parameters users usually need to choose.

For precision, the workflow uses the following internal genotype datasets:

- `raw_sorted.vcf.gz`: the original ipyrad VCF, sorted, compressed, and indexed.
- `subset.vcf.gz`: the raw VCF after optional user-defined sample subsetting.
- `filtered.vcf.gz`: the subset VCF after optional relatedness filtering and variant missing-data filtering.
- `biallelic_snps.vcf.gz`: the filtered VCF restricted to biallelic SNPs with `MAC > 1` and the optional `maf` filter.
- `biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`: the final analysis SNP set after applying `parameters.thinning_strategy`. This is what `get_filtered_vcf_output()` returns and is the main input for most downstream SNP-based analyses.
- `biallelic_snps_thinned.bed/.bim/.fam` and `biallelic_snps_thinned.str`: PLINK and STRUCTURE exports created from the final analysis VCF above. The filenames keep the historical `thinned` basename even when the underlying SNP set actually comes from `ld_pruning` or `parameters.thinning_strategy: none`.

### Population structure and ancestry

**STRUCTURE** uses Bayesian clustering to estimate ancestry proportions and infer the number of genetic clusters (`K`). It uses the STRUCTURE-format export of the final analysis VCF together with `data/mainparams` and `data/extraparams`. `parameters.k_values` is the candidate K range: `K=1` means no subdivision, while larger values allow more clusters. Start broad enough to include plausible structure, inspect model-choice summaries, then narrow if needed. `parameters.structure.replicates` is the number of independent runs per K; more replicates help detect unstable runs and improve consensus Q matrices, but runtime scales roughly linearly.

**fastStructure** is a faster variational-Bayesian alternative to STRUCTURE for ancestry estimation across multiple `K` values. It uses the PLINK export of the final analysis VCF and the same `parameters.k_values` logic as STRUCTURE. `parameters.faststructure.tol` controls convergence tolerance; smaller values ask for stricter convergence and can increase runtime. `parameters.faststructure.prior` can be `"simple"` for a general-purpose prior or `"logistic"` when ancestry is expected to vary with covariates.

**ADMIXTURE** estimates ancestry proportions by maximum likelihood and reports cross-validation error for model choice across `K`. It uses the PLINK export of the final analysis VCF. The main parameter is `parameters.k_values`; include a wide enough range to bracket the likely number of clusters, then use the CV-error and downstream biological interpretation to choose which K values to present.

**sNMF** estimates ancestry proportions by sparse non-negative matrix factorisation (Frichot et al. 2014), implemented in the LEA R package, and is included as a fast, IBD-robust complement to STRUCTURE/ADMIXTURE. It uses the final analysis VCF and `indpopdata.txt`; genotypes are converted to an LEA `.geno` file (alternate-allele dosage `0..ploidy`, missing genotypes coded `9` and handled natively). Because sNMF fits ancestry and ancestral allele frequencies by regularised least squares rather than by an explicit Hardy-Weinberg/linkage-equilibrium likelihood, it is robust to departures from those assumptions — in particular the continuous allele-frequency gradients created by isolation by distance (IBD), which can lead likelihood-based clustering to over-split clinal variation into spurious discrete clusters. It runs once per value in `parameters.k_values`. `parameters.snmf.ploidy` must match your data (`1` for haploids or inbred diploids coded as haploid, `2` for standard diploids). `parameters.snmf.repetitions` sets the number of independent runs per K; sNMF's least-squares optimisation can reach different local optima, so several runs are fit and the one with the lowest cross-entropy is kept (the default of `10` follows the LEA vignette). `parameters.snmf.alpha` is the regularisation parameter on ancestry coefficients: larger values (e.g. `100`) penalise overfitting and are recommended for small or noisy SNP sets, so it is worth comparing a few values. `parameters.snmf.tolerance` and `parameters.snmf.iterations` map to the `snmf()` convergence tolerance and maximum iterations. The number of ancestral populations is chosen from the **cross-entropy criterion** (Alexander & Lange 2011), computed on a masked fraction of held-out genotypes set by `parameters.snmf.percentage` (default `0.05`): the workflow fits all `parameters.k_values`, summarises cross-entropy across repetitions, and produces an Evanno-style cross-entropy-versus-K plot in which the best-supported K is the minimum or the point where cross-entropy plateaus. `parameters.snmf.seed` fixes the random seed for reproducibility, and `parameters.resources.snmf.threads` sets the number of CPUs passed to `snmf()`. Outputs include per-K Q matrices, per-K and summary cross-entropy tables, the cross-entropy choose-K plot, and mapmixture barplots (always) and pie maps (when `popdata` coordinates are available), using the same shared `mapmixture` / `parameters.map_background` styling as the other ancestry methods. Unlike tess3r and conStruct, sNMF is **non-spatial**: pair it with tess3r/conStruct when you want an explicitly spatial view, and use sNMF when you want a fast non-spatial ancestry estimate that remains stable under IBD.

**ALStructure** estimates the same PSD admixture model (ancestry proportions `Q` and ancestral allele frequencies `P`) as STRUCTURE/ADMIXTURE, but fits it with the likelihood-free algorithm of Cabreros & Storey (2019): it first estimates the low-dimensional latent subspace spanned by the ancestry matrix from the spectrum of a bias-corrected genotype covariance matrix, then recovers `P` and `Q` by alternating least squares. Because this is a deterministic spectral procedure rather than MCMC (STRUCTURE) or block-relaxation EM (ADMIXTURE), it has **no random restarts to compare and no local-optima/convergence problems**, and it is comparatively robust to **unbalanced sampling** (very uneven sample sizes across populations bias STRUCTURE-type likelihoods but affect the subspace estimate far less). It is the same model family as SCOPE but installs as a pure-R package (`StoreyLab/alstructure`, fetched once by the `install_alstructure` rule), so it needs no compiler or external binary. It uses the final analysis VCF and `indpopdata.txt`; genotypes are converted to alternate-allele dosages (`0..ploidy`) and missing values are imputed internally by per-locus means. The distinctive feature is **native K estimation**: `estimate_d()` (Leek 2011) estimates the number of ancestral populations from the rank of the latent subspace, and the workflow also writes the latent-subspace eigenvalue scree across `parameters.k_values` and plots it (the elbow corroborates the `estimate_d` value) — so K selection is built into the method rather than requiring a separate cross-validation sweep. `parameters.alstructure.ploidy` must match your data (`1` for haploids/inbred diploids coded as haploid). `parameters.alstructure.svd_method` is `"base"` (exact) or `"truncated_svd"` (faster on very large SNP sets); `parameters.alstructure.tolerance` and `parameters.alstructure.max_iters` control ALS convergence; `parameters.alstructure.order_method` (`"ave_admixture"` or `"var_explained"`) sets the cluster ordering; `parameters.alstructure.seed` fixes the ALS initialisation for reproducibility. ALStructure requires at least two latent populations (`estimate_d` enforces `d >= 2` and the factorisation collapses at `d = 1`), so `K = 1` in `parameters.k_values` is skipped automatically (as for DAPC). Ancestry is fit once per remaining value in `parameters.k_values`, with outputs comprising the `estimate_d` choose-K table and eigenvalue plot, per-K Q matrices, and mapmixture barplots (always) and pie maps (when `popdata` coordinates exist), using the same shared `mapmixture`/`parameters.map_background` styling as the other ancestry methods. Use ALStructure when you want a fast, deterministic, convergence-stable ancestry estimate with an automatic K estimate; it pairs naturally with sNMF (least-squares, IBD-robust) and ADMIXTURE (likelihood, CV-error K) as a methodological triangulation, and with tess3r/conStruct for an explicitly spatial view.

**evalAdmix** evaluates the fit of ancestry models by calculating correlations of residuals after ADMIXTURE, fastStructure, or STRUCTURE inference. It uses the same PLINK export of the final analysis VCF together with the corresponding `Q` and `P` matrices. It has no main numeric parameter in `config.yaml`; use it when you want residual-correlation diagnostics to check whether an inferred admixture model adequately captures covariance among individuals.

**mapmixture** is the shared geographic visualisation layer for ancestry analyses. It uses ancestry coefficients from STRUCTURE, fastStructure, ADMIXTURE, sNMF, ALStructure, tess3r, conStruct, and DAPC together with `indpopdata.txt` to generate population pie-chart maps and barplots. `parameters.mapmixture.pie_size`, `parameters.mapmixture.pie_border`, and `parameters.mapmixture.pie_opacity` control ancestry pies on maps; `parameters.mapmixture.structure_colors` assigns colours to ancestry components, so provide at least as many colours as the largest K you plot. For barplots, `parameters.mapmixture.population_labels` selects one or more `indpopdata.txt` columns used to label/group individuals, `parameters.mapmixture.site_order` sets an explicit population order, and `parameters.mapmixture.population_sort_by` orders populations by a metadata column when `parameters.mapmixture.site_order` is unset. Shared map appearance comes from `parameters.map_background`: `parameters.map_background.basemap` chooses simple/Natural Earth/elevation backgrounds, `parameters.map_background.crs` defines map coordinates, `parameters.map_background.width`/`height`/`dpi` control output size, and `parameters.map_background.arrow`/`scalebar` control cartographic decorations.

**tess3r** estimates spatial ancestry coefficients with geographically constrained non-negative matrix factorisation. It uses the final analysis VCF and `indpopdata.txt`, requires coordinates, and runs once per value in `parameters.k_values`. Genotypes are converted from VCF GT fields to alternate-allele dosage matrices (0..`ploidy`, compatible with the lfmm layout described in the [tess3r vignette](https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html)). `parameters.tess3.ploidy` must match your data (`1` for haploids or inbred diploids coded as haploid, `2` for standard diploids) because `tess3()` has no ploidy default. `parameters.tess3.method` selects the fitting algorithm (`"projected.ls"` is the fast default; `"qp"` is mainly a slower sensitivity check). `parameters.tess3.replicates`, `parameters.tess3.max_iteration`, and `parameters.tess3.tolerance` map to `tess3()` arguments `rep`, `max.iteration`, and `tolerance`. The choose-K run fits all `parameters.k_values` in one call and produces an Evanno-style score plot: by default it summarises training RMSE (`mask: 0`, `crossvalid: false`, `crossentropy: false`, matching `plot.tess3()` defaults); set `parameters.tess3.mask > 0` and `parameters.tess3.crossvalid: true` for masked cross-validation RMSE, or `parameters.tess3.crossentropy: true` to plot cross-entropy instead. Interpolated ancestry surfaces use tess3r's `plot.tess3Q` / `ggtess3Q` defaults unless overridden: `parameters.tess3.map_method` (`"map.max"` or `"map.all"`), `parameters.tess3.interpolation_knots` (`FieldsKrigModel(10)`), and an interpolation grid derived from shared `parameters.map_background` `width`/`height`/`dpi` (rx = width×dpi/10, ry = height×dpi/10). The workflow produces three map styles per K: mapmixture pie maps (shared `parameters.map_background` styling), base-R tess3r interpolation PDFs, and ggplot maps that overlay `ggtess3Q` on the mapmixture basemap (`parameters.tess3.ggmap_point_size` controls sample markers). `parameters.resources.tess3.threads` controls OpenMP cores. Outputs also include per-K Q matrices, cross-entropy tables, tess3r native max-cluster PNGs, and barplots.

**conStruct** models spatial population structure by allowing covariance to decay continuously with geography rather than assuming strictly discrete clusters. It uses the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`) and `indpopdata.txt`, and it requires geographic metadata. It uses `parameters.k_values` for layer counts. `parameters.construct.n_chains` sets the number of independent MCMC chains and `parameters.construct.n_iterations` sets the number of iterations per chain. The conStruct package's own defaults are `1` chain and `1000` iterations; this pipeline deliberately uses larger values — `2` chains and `50000` iterations — because at least two chains are needed to assess convergence (a single chain cannot be checked against another) and more iterations are generally required to describe the posterior surface well, especially with more genotyped samples. These are still only safer starting points: you should examine the model output to assess convergence and increase `parameters.construct.n_iterations` further (and/or run more chains) if convergence is poor. **Warning:** conStruct is computationally heavy and should be run on a small, aggressively thinned SNP set (roughly 1,000-5,000 loci); for RADseq or whole-genome data, thin the markers down specifically for conStruct (and use tess3r on the fuller set for visualisation) or runtimes and memory can become prohibitive. Geographic coordinates and the pairwise distance matrix are taken from `indpopdata` automatically. Layer-proportion maps and barplots use the same shared `mapmixture` / `parameters.map_background` styling as STRUCTURE, ADMIXTURE, tess3r, and DAPC.

Use conStruct when the main question is whether apparent ancestry components remain after accounting for isolation by distance (IBD). This is the key distinction from spatial ancestry methods such as tess3r: tess3r regularizes ancestry coefficients toward spatially smooth solutions, whereas conStruct explicitly models distance-decay of relatedness within each layer and uses extra layers only when the data contain covariance not explained by geography alone. In taxa with strong continuous IBD, ADMIXTURE/STRUCTURE-like approaches can split the ends of a sampling range into different clusters and label central samples as admixed; conStruct is designed to diagnose that artifact by comparing spatial and non-spatial fits, cross-validation scores, and per-layer covariance contributions. In fragmented systems with real barriers, conStruct and spatial ancestry maps should often support similar biological structure, but conStruct gives a more defensible basis for choosing `K`.

conStruct is computationally much heavier than tess3r-style spatial NMF because its Bayesian covariance model evaluates relationships among samples during MCMC. A practical use case is roughly tens to a few hundred individuals and a representative LD-pruned or thinned SNP set, often around 1,000-5,000 loci. For RADseq or whole-genome datasets with many more markers, thin aggressively for conStruct, then use tess3r on the fuller SNP set for visualization, sensitivity checks, or genome scans.

After running the analysis, inspect the per-parameter and posterior-probability trace plots to diagnose MCMC performance: a well-mixed chain looks "fuzzy", whereas a trend that has not plateaued means the chain has not converged and should be run longer (increase `parameters.construct.n_iterations`), and a chain that jumps between two or more levels signals a multi-modal posterior. Also compare the multiple independent chains (`parameters.construct.n_chains`): they should converge on similar parameter estimates and posterior probabilities — if runs disagree but their posterior probabilities are indistinguishable that indicates true multi-modality, while a subset of runs with consistently higher posterior probability indicates the others are mixing poorly. Finally, watch for "homogeneous minimum layer membership" (HMLM), where all samples show suspiciously similar admixture in one layer; this is driven by missing data (e.g. RADseq allelic dropout), so if you see it, reduce missingness by dropping poorly genotyped samples or loci.

**DAPC** (Discriminant Analysis of Principal Components) provides multivariate clustering and assignment based on the final analysis VCF. `parameters.k_values` defines the cluster counts tested, with `K=1` skipped because DAPC requires at least two groups. `parameters.dapc.n_pca` is the number of PCA axes retained before discrimination: too few can miss structure, while too many can overfit noise. `parameters.dapc.n_da` is the number of discriminant axes retained, usually bounded by K minus one. `parameters.dapc.criterion` controls automatic cluster selection; `"diffNgroup"` is a practical default when criterion curves are not obvious.

**sPCA** (Spatial PCA, adegenet 2.1.x) ordinates individuals using axes that combine genetic variance with spatial autocorrelation (Moran's eigenvector maps). It requires spatial metadata in `indpopdata.txt` and follows the [adegenet sPCA tutorial](https://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf) for global/local tests, score plots, and loading plots. `parameters.spca.n_pca` optionally reduces the genotype matrix before sPCA (`null` uses all loci). `parameters.spca.nfposi` and `parameters.spca.nfnega` set retained global and local axes: positive/global axes describe broad clines, while negative/local axes describe fine-scale patches; tune them from the eigenvalue plots. `parameters.spca.type` chooses the neighbour graph: type `1` Delaunay is a flexible default, type `6` K-nearest neighbours needs `k`, and type `5` distance bands need `d1`/`d2`. `parameters.spca.nperm` controls spatial tests; with `999`, the smallest possible p-value is about `0.001`. `parameters.spca_plot.color_by`, `parameters.spca_plot.pc_max`, and `parameters.spca_plot.loading_threshold` control score/loadings plots, while `parameters.spca_map` controls point size, opacity, and score colours on repo-style maps.

**OrientAGraph** infers a population graph with optional migration edges from population allele counts using TreeMix-compatible inputs and outputs. It uses the final analysis VCF and `indpopdata.txt`; `parameters.treemix.population_column` selects the metadata column used to group individuals into populations. The converter writes the standard TreeMix `ref,alt` allele-count matrix per population, along with a `.clust` file, retained SNP positions, and a population label mapping. `parameters.treemix.migration_edges` lists the `-m` values to run, `parameters.treemix.k` controls the SNP block size (`-k`), and `parameters.treemix.root`, `parameters.treemix.global`, and `parameters.treemix.seed` are passed through to OrientAGraph when set. Set `parameters.treemix.root` to a real outgroup/population label from `treemix.populations.tsv` for final analyses; leaving `parameters.treemix.root` unset (`null`) is mainly for exploratory unrooted runs. Maximum likelihood network orientation is enabled by default with `treemix.mlno: true`, which passes bare `-mlno` for the pinned bioconda build; set `false` or `null` to omit it. If using an OrientAGraph 1.2+ source build, `parameters.treemix.mlno` can instead be a value such as `"1,2"`, `"0"`, or `"all"`. `parameters.treemix.allmigs` controls exhaustive legal migration-edge search and defaults to `false` because the original-paper mode can be very slow. `parameters.treemix.bootstrap` is an opt-in support workflow, not part of default graph fitting: set `parameters.treemix.bootstrap.enabled: true`, choose final `parameters.treemix.bootstrap.migration_edges` such as `[3]`, and use `parameters.treemix.bootstrap.replicates: 100` for final support checks. Each bootstrap job passes `-bootstrap` once with an independent deterministic seed, then the workflow summarizes final likelihoods and migration-edge support across replicates. Plot outputs include original TreeMix `plot_tree()` graph PDFs, original `plot_resid()` residual PDFs, an OptM migration-edge selection table/plot, a simple likelihood-by-migration-count summary plot, and optional bootstrap likelihood/support summaries. `parameters.treemix.optm.method` selects the OptM method (`Evanno`, `linear`, or `SiZer`), `parameters.treemix.optm.replicates` controls independent runs per `m` for OptM, and `parameters.treemix.plot` controls plot dimensions and resolution. Following Fitak (2021), OptM is most informative with multiple independent runs per migration-edge count (for example `10+` for final analyses) and more than three usable migration-edge levels; if higher requested `m` runs cannot add additional edges, the workflow falls back to a likelihood summary and reports requested versus achieved migration counts.

### Ordination and genetic distances

**PCAone** performs principal component analysis on the PLINK export of the final analysis VCF and is intended as the main fast PCA implementation for large SNP matrices. `parameters.PCAone.PCnum` controls how many PCs are calculated, while `parameters.pca_plot.pc_max` controls how many axes are plotted in pairwise combinations. `parameters.PCAone.SVD_method` selects the PCA algorithm: `3` is full SVD and accurate for smaller datasets, `2` is PCAone's randomized window-based method for large datasets, and `0`/`1` are faster approximate methods. `parameters.pca_plot.color_by` lists metadata columns used to colour/facet plots, and `parameters.pca_plot.colors` provides a manual categorical palette.

**PCAone with EMU** uses the EMU mode within PCAone to model missing genotypes during PCA. It uses the same inputs and output structure as standard PCAone, but is preferable when missingness is substantial and should be handled during ordination rather than by prior filtering. It uses the same `parameters.PCAone.PCnum` and `parameters.pca_plot` settings as PCAone. `parameters.PCAone.EMU_SVD_method` selects the SVD method for EMU mode (`0`–`2`; do not use full SVD `3` here). The default is **`2`** (window-based randomized SVD / winSVD), which is stable and fast for typical RAD datasets; `0` (IRAM in-core) can segfault on some inputs.

**PCAone with missing-data thresholds** repeats PCAone after creating separate SNP datasets filtered at each value in `parameters.PCAone.miss`. This allows explicit sensitivity analysis of ordination results to marker missingness. Each threshold is run independently, so long `miss` lists multiply runtime and result sets.

**EMU-PCA** runs the standalone `emu` program rather than PCAone's internal EMU mode. It is intended for low-coverage or uncertainty-prone genotype data where expectation-maximisation-based PCA is preferred. It uses the PLINK export of the final analysis VCF and the same plotting family as PCAone. `parameters.EMU.eig` selects the first PC reported by EMU, `parameters.EMU.eig_out` controls how many PCs are written, and `parameters.pca_plot` controls colours and plotted axis combinations.

**VCF2PCACluster** performs PCA with its own internal SNP filtering and optional clustering/kinship correction, using the final analysis VCF as input. It starts from the main downstream SNP set and then applies its own internal `MAF`, missingness, heterozygosity, and HWE filters. `parameters.vcf2pcacluster.PCnum` sets the number of PCs, `parameters.vcf2pcacluster.cluster_method` controls optional clustering (`Kmean`, `EM`, `DBSCAN`, or `None`), and `parameters.vcf2pcacluster.SNP_filtering` is separate from the workflow's main preprocessing filters, making it useful for method-specific sensitivity analyses. `parameters.vcf2pcacluster.KinshipMethod` selects the tool's kinship correction estimator; keep the default unless you know which estimator fits your design.

**Euclidean genetic distance** (and the Kosman-Leonard and average-squared distances) is computed from the PLINK export of the **un-thinned biallelic SNP set** (`biallelic_snps.vcf.gz`), treating genotypes as alternate-allele dosages (`0`, `1`, `2`) and mean-imputing missing values per SNP. Distance-based methods do not assume unlinked markers, so the full biallelic SNP set is used instead of the thinned set — this retains more data and makes the result deterministic (the thinned set depends on random tie-breaking). The set still carries the `MAC > 1` (and optional `maf`) filtering from `select_biallelic_snps`. This matrix is the downstream input for PCoA and MAPI.

**PCoA** performs principal coordinates analysis on the Euclidean distance matrix and reuses the same plotting system as PCAone for colored, labeled, missingness-based, and faceted ordinations. `parameters.pca_plot.pc_max` controls plotted axis combinations, `parameters.pca_plot.color_by` selects metadata columns for colouring/faceting, and `parameters.pca_plot.colors` supplies categorical colours.


### Genetic diversity, relatedness, and demographic summaries

**Pixy** calculates unbiased estimates of nucleotide diversity (`pi`), genetic differentiation (`FST`), and sequence divergence (`dXY`) using both variant and invariant sites extracted from the original ipyrad `.loci` file rather than from the downstream SNP-only VCFs. It uses the project keep-list after sample filtering plus grouping information from `indpopdata.txt`. `parameters.pixy.stats` chooses which statistics to compute: `pi` is within-population diversity, `fst` is relative differentiation, and `dxy` is absolute divergence. `parameters.pixy.window_size` sets genomic window size; smaller windows give finer but noisier resolution, while larger windows smooth estimates. `parameters.pixy.bootstrap_replicates` controls uncertainty estimates. `parameters.pixy.group_by` defines the metadata columns used as pixy populations; include `Site` if you want per-site pi summaries and geographic pi maps. `parameters.pixy.arrange_by` changes grouped/sorted pi barplot variants but not the pixy calculation itself. Summary files are constructed from raw windows using weighted means (`no_sites` for `pi`/`dXY`, `no_snps` for `FST`) and bootstrap confidence intervals.

**Genome scan** performs sliding-window analyses of `FST`, `pi`, and `dXY` between two user-defined groups. It uses an invariant-site VCF reconstructed from the original ipyrad data, together with `indpopdata.txt`, to identify samples whose `parameters.genome_scan.population_column` value matches `parameters.genome_scan.pop1` or `parameters.genome_scan.pop2`; these values must exactly match metadata entries. The workflow subsets the invariant-site VCF to those two groups, creates a two-population popmap, and runs `pixy` on the reduced dataset, so samples outside the focal comparison are excluded. `parameters.genome_scan.window_size_fst`, `parameters.genome_scan.window_size_pi`, and `parameters.genome_scan.window_size_dxy` set statistic-specific window sizes; FST often uses smaller windows than pi/dXY in the default config. Note that for de novo RAD data each locus is a short, independent contig with no inter-locus ordering, so these "sliding windows" are effectively per-RAD-locus statistics rather than a true genomic scan; genomic-position interpretation (e.g. clustering of differentiation peaks) only applies when loci are mapped to a reference genome.

**Relatedness** estimates pairwise kinship using multiple estimators, including the Yang et al. genomic relatedness statistic, Manichaikul/KING-style kinship, PLINK IBD summaries, and PLINK2 KING outputs. The main relatedness analysis uses the final analysis VCF and its derived PLINK files, whereas the earlier relatedness-filtering step in preprocessing uses `subset.vcf.gz` before variant missingness filtering to decide which individuals to retain. `parameters.relatedness_plot.color_by` selects metadata columns used to colour nodes (`"none"` disables colouring), `parameters.relatedness_plot.plot_all: false` plots only related individuals, and `parameters.relatedness_plot.plot_all: true` also includes unrelated individuals as isolated nodes. `parameters.relatedness_plot.colors` can provide a manual palette.

**ROH** uses `bcftools roh` to identify runs of homozygosity and summarises them across individuals and groups. It uses `filtered.vcf.gz` (the variant set after sample subsetting, relatedness filtering, and the optional `parameters.vcf_filtering.f_missing` missingness filter, before the biallelic/`MAC > 1`/`maf` and thinning filters) together with `indpopdata.txt` — the **same filtered dataset used by fineRADstructure**. This un-thinned set is used deliberately: ROH detection relies on the density of consecutive genotypes along each RAD locus, so thinning to one SNP per locus would remove the consecutive sites that ROH calling requires. `parameters.roh.group_by` lists metadata columns used to summarise and compare ROH patterns. ROH segments are classified by length as long (>5 Mb; ≤~10 generations), medium (1-5 Mb; ~10-50 generations), or short (<1 Mb; >~50 generations), using G = 100/(2 × cM) with centimorgans estimated from megabases assuming 1 cM/Mb. When `popdata` is configured with coordinates, the workflow also builds a geographic map of mean per-site `F_ROH` using the same map styling as pixy π maps (`parameters.roh.point_size`, `parameters.roh.map_outline`, and shared `parameters.map_background` settings). Note that ROH calling has limited power on de novo RAD data, where each locus is a short, independent contig with no inter-locus ordering; it is most meaningful when loci are mapped to a reference genome so that long runs of consecutive homozygous sites can actually be detected.

**AMOVA** partitions molecular variance across hierarchical levels defined in the metadata. It uses the final analysis VCF (`biallelic_snps_{thinned|ld_pruned|all_snps}.vcf.gz`) and `indpopdata.txt`. `parameters.amova.strata` should list metadata levels from broad to fine, for example `["Region", "Site"]`; `parameters.amova.nperm` controls permutation tests, with more permutations giving more precise p-values.

### Spatial landscape genetics

**MAPI** (Mapping Averaged Pairwise Information) identifies spatial regions where observed genetic distances are unusually high or low relative to geography. It uses the Euclidean distance matrix and `indpopdata.txt`. `parameters.mapi.n_permutations` controls permutation testing, `parameters.mapi.grid_halfwidth` sets hexagonal grid-cell half-width in metres, and smaller cells give finer but noisier/slower maps. `parameters.mapi.crs_projected` should be suitable for distance calculations in the study area, `parameters.mapi.crs_geographic` describes the input coordinate CRS, and `parameters.mapi.alpha` is the significance cutoff for upper/lower tails. `parameters.mapi.fill_var`, point/tail styling, and tail colours control the final map appearance.

### Phylogenetic and historical inference

**NeighborNet (fastnntr)** infers a split network from the p-distance matrix using `fastnntr::run_neighbournet_networkx` (Rust-backed NeighborNet; phangorn is retained for planar coordinates and the reference `plot.networx` figure). The p-distance is computed from `filtered.vcf.gz` — the missingness-filtered variant set after sample subsetting and optional relatedness filtering, **before** biallelic/`MAC > 1`/`maf` and thinning filters (the same input stage as fineRADstructure). It is the pairwise mean of `|g_i - g_j| / 2` over biallelic variant sites genotyped in both individuals (multiallelic sites skipped). Networks are drawn with `tanggle::ggsplitnet()` (the standard phangorn/ggtree split-network renderer), with tip colours from metadata and `geom_tiplab2()` labels; phangorn's reference `plot.networx` figure is also written. `parameters.neighbornet.ordering_method` and `parameters.neighbornet.inference_method` select the fastnntr cycle-ordering and split-weight solvers (defaults: `closest-pair` + `active-set`, matching SplitsTree6-style settings). `parameters.neighbornet.color_by` selects metadata columns used to colour tips; per-column tip colours follow the same `parameters.neighbornet.<column>.colors` pattern as ROH group plots. `parameters.neighbornet.width`, `parameters.neighbornet.height`, and `parameters.neighbornet.dpi` set plot dimensions/resolution. `parameters.neighbornet.linewidth` sets split-network edge width for both the ggsplitnet and phangorn figures (default `0.5`).

**IQ-TREE** performs maximum-likelihood phylogenetic inference from the ipyrad `.phy` alignment. `parameters.iqtree.model` sets the substitution model; `"MFP"` asks IQ-TREE ModelFinder Plus to choose one. `parameters.iqtree.bootstraps` controls ultrafast bootstrap replicates, where more replicates improve support stability but increase runtime. Rooting uses `parameters.iqtree.outgroup` when provided (single sample or comma-separated sample set), otherwise midpoint rooting is used automatically. `parameters.iqtree.robust-phy` controls the workflow's robust-tree selection, and `parameters.iqtree.support_threshold` is the minimum bootstrap support displayed on plotted trees.

**trimAl + IQ-TREE** adds an alignment-trimming step before phylogenetic inference. It uses `trimAl` to remove sites with excessive gaps and then runs IQ-TREE on the trimmed alignment. `parameters.trimal.gt` is the gap threshold: values close to `1` are permissive, while lower values remove columns with more gaps. The downstream tree step uses the same `iqtree` settings.

**njtree** runs [rapidNJ](https://github.com/somme89/rapidNJ) on the ipyrad concatenated alignment. It supports nucleotide distance models (`parameters.njtree.distance_model`: `kim` for Kimura 2-parameter, `jc` for Jukes-Cantor) and optional bootstrap replicates (`parameters.njtree.bootstraps`; set to `0` to skip). This requires `{ipyrad_prefix}.phy` to exist in the ipyrad output directory — if your ipyrad run did not export PHYLIP, re-export from ipyrad first. The workflow converts `.phy` to FASTA internally (rapidNJ cannot read PHYLIP alignments directly). `parameters.iqtree.support_threshold` controls which bootstrap values are shown on the unrooted ggtree plot.

**fineRADstructure** uses `RADpainter` and `fineRADstructure` to infer recent coancestry and hierarchical relationships among individuals from RAD-seq loci. It takes `filtered.vcf.gz`, that is, the missingness-filtered variant set after sample subsetting and relatedness filtering, together with `indpopdata.txt`. This set deliberately keeps **all** SNPs — including singletons and multiallelic sites — and is **not** passed through the `MAC > 1`, `maf`, or one-SNP-per-locus thinning filters used by the other analyses. fineRADstructure's recent-coancestry signal is dominated by rare alleles, so removing singletons (as `biallelic_snps.vcf.gz` does) would discard the most informative variants; the un-thinned set also keeps multiple SNPs per RAD locus so haplotypes can be reconstructed. The only variant filter applied is the optional `parameters.vcf_filtering.f_missing` missingness threshold (disabled by default); this is a site-level filter that drops poorly-genotyped loci, not rare alleles, so it does not undercut the rare-allele signal described above. `parameters.fineradstructure.cluster.mcmc_iterations` controls clustering MCMC length, `parameters.fineradstructure.cluster.burnin` is discarded before summarising the chain, and `parameters.fineradstructure.cluster.thinning` controls how often samples are retained. `parameters.fineradstructure.tree.mcmc_iterations` controls tree-building MCMC. `parameters.fineradstructure.plot.max_indv` and `parameters.fineradstructure.plot.max_pop` cap heatmap size for large datasets. Final coancestry heatmaps are written under `results/<project>/fineradstructure/plots/` (`SimpleCoancestry.pdf`, `PopAveragedCoancestry.pdf`, `PopAveragedCoancestryLabeled.pdf`, plus a `.fineradstructure_plots.rds` object for replotting).

### Population metadata

**Population assignment** is generated automatically from either a user-provided `popmap`, a richer `popdata` table with geographic coordinates and grouping variables, or a `popseparator` that parses population codes from sample names. The main inputs are `parameters.popmap`, `parameters.popdata`, and `parameters.popseparator`, and the main outputs are `results/<project>/indpopdata.txt`, `results/<project>/indpopdata_all.txt`, and `results/<project>/stats_samples/<project>.population_summary.txt`.

**Population maps and summaries** use `indpopdata.txt` and the shared map settings to create descriptive sampling maps. If `popdata` is available, the pipeline can generate `results/<project>/stats_samples/plots/<project>.population_map.pdf` together with tabular summaries of sample counts per population. `parameters.population_map` controls point size/shape/colour and label display (`parameters.population_map.show_points`, `parameters.population_map.show_labels`, `parameters.population_map.label_size`, `parameters.population_map.label_color`, `parameters.population_map.label_fontface`), ggrepel label placement (`parameters.population_map.force`, `parameters.population_map.force_pull`, `parameters.population_map.min.segment.length`, `parameters.population_map.segment.color`, `parameters.population_map.segment.size`), and `parameters.population_map.max.overlaps`, which sets how many overlapping site labels ggrepel keeps: the default `".inf"` shows all labels, while smaller integers suppress crowded labels with more than that many overlaps. Keep `max.overlaps` quoted in YAML — unquoted `.inf` is parsed as floating-point infinity and breaks Snakemake's R script preambles. `parameters.map_boundary` can crop all maps to the same extent, and `parameters.map_background` controls basemap type, CRS, output size, resolution, north arrow, and scale bar.

### Computational resources

Resources per analysis type defined under `parameters.resources`:
- `parameters.resources.<module>.threads` <INTEGER>: Number of CPU cores
- `parameters.resources.<module>.mem_mb` <INTEGER>: Memory allocation in megabytes
- `parameters.resources.<module>.runtime` <INTEGER>: Maximum runtime in minutes