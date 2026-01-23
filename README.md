# phylogeo-snakemake

Snakemake pipeline for population genetics and phylogeographic data analysis from ipyrad output (`.vcf` and `.loci` files).

## Configuration

All settings are in `config/config.yaml`:
- Define multiple projects under `config["projects"]`, each with its own `ipyrad_prefix` and analysis selection
- Global parameters (under top-level `parameters:`) are inherited by all projects unless overridden in `config["projects"][project_name]["parameters"]`
- Enable/disable analyses per project in `config["projects"][project_name]["analyses"]`

## Running

Best is to use conda/mamba:
```
snakemake --use-conda --conda-frontend mamba
```

On HPC you need to specify the executor and resources, eg.:

```
snakemake --executor slurm --default-resource slurm_account=plgdryas-cpu slurm_partition=plgrid  --use-conda --conda-frontend mamba --jobs 10
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
- **Parameters**: None (uses sample list from config)

#### 3. Relatedness Filtering (Optional)
- **KING-based filtering**: Filter related/clonal individuals using plink2 KING method
- **Parameters** (from `config["projects"][project_name]["parameters"]["relatedness_filtering"]`):
  - `enabled` <BOOLEAN>: Enable/disable relatedness filtering (default: `false`)
  - `king_threshold` <FLOAT>: KING kinship threshold for filtering (default: `0.0884`, equivalent to PI-HAT > 0.2, 2nd-degree relatives)
    - Common values: `0.354` (duplicates/twins), `0.177` (1st-degree), `0.0884` (2nd-degree), `0.0442` (3rd-degree)
  - `mac_threshold` <INTEGER>: Minimum Minor Allele Count (MAC) for KING calculation (default: `1`)
  - `group_by` <LIST[STRING]>: List of columns from popdata to group barplots by (default: `["Region"]`)
- **Output**: Generates filtered sample list and KING relatedness table; creates barplots showing removed individuals grouped by specified columns

#### 4. Missing Data Filtering (After Relatedness)
- **Variant-level missing data filter**: Applied after relatedness filtering (if enabled)
- **Parameters** (from `config["projects"][project_name]["parameters"]["vcf_filtering"]`):
  - `f_missing` <FLOAT>: Maximum missing data threshold for variant filtering (default: `0.5`)
    - Variants missing in > `f_missing` proportion of samples are filtered out
    - `1` keeps all variants, `0` removes variants with any missing data
- **Note**: If `f_missing >= 1`, missing data filtering is skipped

#### 5. Biallelic SNP Filtering
- **Hard-coded filter**: Retains only biallelic SNPs with Minor Allele Count (MAC) > 1
- **Parameters**: None (hard-coded: `MAC > 1`, `-m2 -M2`, `-v snps`)

#### 6. VCF Thinning (One SNP per RAD Fragment)
- **RAD fragment thinning**: Select one SNP per RAD fragment to ensure unlinked markers
- **Parameters** (from `config["projects"][project_name]["parameters"]["vcf_thinning"]`):
  - `min_coverage` <INTEGER>: Minimum number of samples with data at a locus (default: `0`)
  - `method` <STRING>: SNP selection method - `random`, `max_coverage`, or `weighted` (default: `max_coverage`)
  - `ties` <STRING>: How to handle tied coverage values - `first` or `random` (default: `random`)
  - `ns_tag` <STRING>: VCF INFO field tag for number of samples with data (default: `NS`)
  - `id_pattern` <STRING>: Regex pattern to extract RAD fragment ID from variant names (default: `loc(\d+)_`)

#### 7. Format Conversion
- **PLINK format**: Export thinned VCF to PLINK format (`.bed/.bim/.fam`)
- **Structure format**: Export thinned VCF to Structure format (`.str`)
- **Missing data filtered PLINK files**: Generate PLINK files filtered by missing data thresholds defined in `config["projects"][project_name]["parameters"]["PCAone"]["miss"]` for PCAone robustness analysis

#### Summary of Filtering Order
1. Sort and index input VCF
2. (Optional) Subset to user-specified samples
3. (Optional) Filter related individuals using KING method
4. Filter variants by missing data threshold (`f_missing`)
5. Filter to biallelic SNPs with MAC > 1
6. Thin to one SNP per RAD fragment
7. Export to various formats (PLINK, Structure)
8. (For PCAone) Generate additional VCFs filtered by multiple missing data thresholds

### Population structure

**Structure** (enabled with `config["projects"][project_name]["analyses"]["structure"]: true`)

Bayesian clustering analysis using STRUCTURE. Requires `data/mainparams` and `data/extraparams` configuration files - run parameters can be tweaked by modifying these files. Numbers of individuals and loci are automatically inferred from .ustr file so this does not need to be changed.

Parameters (from `config["projects"][project_name]["parameters"]`):
- `k_values` <LIST[INTEGER]>: List of K values (number of populations) to test (e.g., `[1, 2, 3]`)

Parameters (from `config["projects"][project_name]["parameters"]["structure"]`):
- `replicates` <INTEGER>: Number of independent runs per K value (default: 2)

**fastStructure** (enabled with `config["projects"][project_name]["analyses"]["faststructure"]: true`)

Fast variational inference alternative to STRUCTURE.

Parameters (from `config["projects"][project_name]["parameters"]`):
- `k_values` <LIST[INTEGER]>: List of K values to test

Parameters (from `config["projects"][project_name]["parameters"]["faststructure"]`):
- `tol` <FLOAT>: Convergence criterion (default: `10e-6`)
- `prior` <STRING>: Prior model - `simple` or `logistic` (default: `simple`)

**Admixture** (enabled with `config["projects"][project_name]["analyses"]["admixture"]: true`)

Maximum likelihood estimation of ancestry proportions. Outputs cross-validation errors for optimal K selection.

Parameters (from `config["projects"][project_name]["parameters"]`):
- `k_values` <LIST[INTEGER]>: List of K values to test

### Principal Component Analysis

**PCAone** (enabled with `config["projects"][project_name]["analyses"]["pcaone"]: true`)

Fast PCA on complete dataset without missing data filtering.

Parameters (from `config["projects"][project_name]["parameters"]["PCAone"]`):
- `SVD_method` <INTEGER>: Singular value decomposition method - 0=IRAM, 1=randomized, 2=window-based, 3=full SVD (default: 3)
- `PCnum` <INTEGER>: Number of principal components to compute (default: 10)

**PCAone with EMU** (enabled with `config["projects"][project_name]["analyses"]["pcaone_emu"]: true`)

PCA using the EMU (Expectation-Maximization-with-Uncertainty) method for handling missing data. Uses same parameters as PCAone above.

**PCAone with missing filters** (enabled with `config["projects"][project_name]["analyses"]["pcaone_miss"]: true`)

Run PCA on datasets filtered by different missing data thresholds to assess impact of data completeness. **Note**: This creates separate VCF files filtered by each missing data threshold in the `miss` list, then runs PCAone on each filtered dataset.

Parameters (from `config["projects"][project_name]["parameters"]["PCAone"]`):
- `miss` <LIST[FLOAT]>: List of maximum missing data thresholds to test (e.g., `[0.1, 0.2, 0.3]`)
  - Each threshold filters variants where `F_MISSING < threshold`
  - Separate PCAone runs are performed for each threshold
- `SVD_method`, `PCnum`: Same as PCAone

**VCF2PCACluster** (enabled with `config["projects"][project_name]["analyses"]["vcf2pcacluster"]: true`)

PCA with kinship correction and flexible SNP filtering. **Note**: VCF2PCACluster applies its own SNP filtering parameters independently from the main preprocessing pipeline. It uses the thinned VCF as input but applies additional filters.

Parameters (from `config["projects"][project_name]["parameters"]["vcf2pcacluster"]`):
- `cluster_method` <STRING>: Clustering algorithm - `EM`, `Kmean`, `DBSCAN`, or `None` (default: `Kmean`)
- `PCnum` <INTEGER>: Number of principal components (default: 10)
- `KinshipMethod` <INTEGER>: Kinship calculation method - 1=Normalized_IBS (Yang/BaldingNicols), 2=Centered_IBS (VanRaden), 3=IBSKinshipImpute, 4=IBSKinship, 5=p_dis (default: 1)

**VCF2PCACluster SNP Filtering** (from `config["projects"][project_name]["parameters"]["vcf2pcacluster"]["SNP_filtering"]`):
- `MAF` <LIST[FLOAT]>: List of minimum minor allele frequency thresholds - multiple thresholds tested (e.g., `[0, 0.001, 0.01, 0.05]`)
- `Miss` <LIST[FLOAT]>: List of maximum missing data thresholds - multiple thresholds tested (e.g., `[0.3, 0.4, 0.5]`)
- `Het` <FLOAT>: Maximum heterozygosity ratio filter (default: `1.0`)
- `HWE` <FLOAT>: Hardy-Weinberg Equilibrium exact test p-value threshold (default: `0`, no filtering)
- `Fchr` <STRING>: Filter specific chromosome (default: empty string for no filter)
- **Note**: All combinations of `MAF` and `Miss` values are tested, generating separate PCA results for each combination

**PCA plotting**

Parameters (from `config["projects"][project_name]["parameters"]["pca_plot"]`):
- `pc1` <STRING>: First principal component axis to plot (default: `"1"`)
- `pc2` <STRING>: Second principal component axis to plot (default: `"2"`)
- `color_by` <LIST[STRING]>: List of popdata columns to use for coloring points (e.g., `["Population", "Region"]`)

### Genetic diversity

**Pixy** (enabled with `config["projects"][project_name]["analyses"]["pixy"]: true`)

Calculate unbiased nucleotide diversity (π), FST, and dXY using both variant and invariant sites extracted from the `.loci` file.

Parameters (from `config["projects"][project_name]["parameters"]["pixy"]`):
- `stats` <STRING>: Space-separated statistics to compute - options: `pi`, `fst`, `dxy` (default: `"pi fst dxy"`)
- `window_size` <INTEGER>: Sliding window size in base pairs (default: 10000)

### Relatedness

**Relatedness** (enabled with `config["projects"][project_name]["analyses"]["relatedness"]: true`)

Calculate pairwise relatedness between individuals using two methods:
- Yang et al. (2010) method: Ajk statistic (expectation: 1 for self, 0 for unrelated within population)
- Manichaikul et al. (2010) method: Kinship coefficient (expectation: 1/2 for monozygotic twins, 1/4 for full-sibs/parent-offspring, 1/8 for second degree, 1/16 for third degree, 0 for unrelated)

### Population metadata

**Population assignment**

Parameters (from `config["projects"][project_name]["parameters"]`):
- `popdata` <STRING>: Path to tab-delimited file with columns: population, latitude, longitude, and optional grouping variables (e.g., region)
- `popseparator` <STRING>: Alternative to popdata - character separating individual ID from population code in sample names (e.g., `"-"` for `ind1-pop1`)

### Computational resources

Resources per analysis type defined under `config["projects"][project_name]["parameters"]["resources"]`:
- `threads` <INTEGER>: Number of CPU cores
- `mem_mb` <INTEGER>: Memory allocation in megabytes
- `runtime` <INTEGER>: Maximum runtime in minutes

## Filtering Parameters Summary

This section provides a quick reference for all filtering parameters used throughout the pipeline.

### Main Preprocessing Filters

| Filter Step | Parameter Path | Parameter | Type | Default | Description |
|------------|----------------|-----------|------|---------|-------------|
| Sample subsetting | `config["projects"][project]["samples"]` | `samples` | STRING (comma-separated) | - | Optional: comma-separated list of sample IDs to retain |
| Relatedness filtering | `relatedness_filtering.enabled` | `enabled` | BOOLEAN | `false` | Enable/disable KING-based relatedness filtering |
| Relatedness filtering | `relatedness_filtering.king_threshold` | `king_threshold` | FLOAT | `0.0884` | KING kinship threshold (0.354=twins, 0.177=1st-degree, 0.0884=2nd-degree, 0.0442=3rd-degree) |
| Relatedness filtering | `relatedness_filtering.mac_threshold` | `mac_threshold` | INTEGER | `1` | Minimum MAC for KING calculation |
| Missing data | `vcf_filtering.f_missing` | `f_missing` | FLOAT | `0.5` | Maximum missing data threshold (1=keep all, 0=no missing allowed) |
| VCF thinning | `vcf_thinning.min_coverage` | `min_coverage` | INTEGER | `0` | Minimum number of samples with data at a locus |
| VCF thinning | `vcf_thinning.method` | `method` | STRING | `"max_coverage"` | SNP selection: `random`, `max_coverage`, or `weighted` |
| VCF thinning | `vcf_thinning.ties` | `ties` | STRING | `"random"` | Handle ties: `first` or `random` |
| VCF thinning | `vcf_thinning.ns_tag` | `ns_tag` | STRING | `"NS"` | VCF INFO field tag for number of samples |
| VCF thinning | `vcf_thinning.id_pattern` | `id_pattern` | STRING | `"loc(\\d+)_"` | Regex pattern to extract RAD fragment ID |

### Analysis-Specific Filters

| Analysis | Parameter Path | Parameter | Type | Default | Description |
|----------|----------------|-----------|------|---------|-------------|
| PCAone missing | `PCAone.miss` | `miss` | LIST[FLOAT] | `[0.1, 0.2, 0.3, ...]` | List of missing data thresholds for robustness testing |
| VCF2PCACluster | `vcf2pcacluster.SNP_filtering.MAF` | `MAF` | LIST[FLOAT] | `[0, 0.001, 0.01, 0.05]` | List of minimum minor allele frequency thresholds |
| VCF2PCACluster | `vcf2pcacluster.SNP_filtering.Miss` | `Miss` | LIST[FLOAT] | `[0.3, 0.4, 0.5]` | List of maximum missing data thresholds |
| VCF2PCACluster | `vcf2pcacluster.SNP_filtering.Het` | `Het` | FLOAT | `1.0` | Maximum heterozygosity ratio filter |
| VCF2PCACluster | `vcf2pcacluster.SNP_filtering.HWE` | `HWE` | FLOAT | `0` | Hardy-Weinberg Equilibrium p-value threshold (0=no filter) |
| VCF2PCACluster | `vcf2pcacluster.SNP_filtering.Fchr` | `Fchr` | STRING | `""` | Filter specific chromosome (empty=no filter) |

