# phylogeo-snakemake

Snakemake pipeline for population genetics and phylogeographic data analysis from ipyrad output (`.vcf` and `.loci` files).

## Configuration

All settings are in `config/config.yaml`:
- Define multiple projects under `config["projects"]`, each with its own `ipyrad_prefix` and analysis selection
- Global parameters (under top-level `parameters:`) are inherited by all projects unless overridden in `config["projects"][project_name]["parameters"]`
- Enable/disable analyses per project in `config["projects"][project_name]["analyses"]`

## Workflow

### Preprocessing

**Sort and filter VCF** - Sort, compress, and index ipyrad VCF output, then filter for biallelic SNPs with MAC>1.

**Thin to one SNP per RAD fragment**

Parameters (from `config["projects"][project_name]["parameters"]["vcf_thinning"]`):
- `min_coverage` <INTEGER>: Minimum number of samples with data at a locus (default: 0)
- `method` <STRING>: SNP selection method - `random`, `max_coverage`, or `weighted` (default: `max_coverage`)
- `ties` <STRING>: How to handle tied coverage values - `first` or `random` (default: `random`)
- `ns_tag` <STRING>: VCF INFO field tag for number of samples with data (default: `NS`)
- `id_pattern` <STRING>: Regex pattern to extract RAD fragment ID from variant names (default: `loc(\d+)_`)

**Format conversion** - Export thinned VCF to PLINK format (`.bed/.bim/.fam`) and Structure format (`.str`). Generate PLINK files filtered by missing data thresholds defined in `config["projects"][project_name]["parameters"]["PCAone"]["miss"]`.

### Population structure

**Structure** (enabled with `config["projects"][project_name]["analyses"]["structure"]: true`)

Bayesian clustering analysis using STRUCTURE. Requires `data/mainparams` and `data/extraparams` configuration files.

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

Run PCA on datasets filtered by different missing data thresholds to assess impact of data completeness.

Parameters (from `config["projects"][project_name]["parameters"]["PCAone"]`):
- `miss` <LIST[FLOAT]>: List of maximum missing data thresholds to test (e.g., `[0.1, 0.2, 0.3]`)
- `SVD_method`, `PCnum`: Same as PCAone

**VCF2PCACluster** (enabled with `config["projects"][project_name]["analyses"]["vcf2pcacluster"]: true`)

PCA with kinship correction and flexible SNP filtering.

Parameters (from `config["projects"][project_name]["parameters"]["vcf2pcacluster"]`):
- `cluster_method` <STRING>: Clustering algorithm - `EM`, `Kmean`, `DBSCAN`, or `None` (default: `Kmean`)
- `PCnum` <INTEGER>: Number of principal components (default: 10)
- `KinshipMethod` <INTEGER>: Kinship calculation method - 1=Normalized_IBS (Yang/BaldingNicols), 2=Centered_IBS (VanRaden), 3=IBSKinshipImpute, 4=IBSKinship, 5=p_dis (default: 1)

Parameters (from `config["projects"][project_name]["parameters"]["vcf2pcacluster"]["SNP_filtering"]`):
- `MAF` <LIST[FLOAT]>: List of minimum minor allele frequency thresholds (e.g., `[0, 0.001, 0.01, 0.05]`)
- `Miss` <LIST[FLOAT]>: List of maximum missing data thresholds (e.g., `[0.3, 0.4, 0.5]`)
- `Het` <FLOAT>: Maximum heterozygosity ratio filter (default: 1.0)
- `HWE` <FLOAT>: Hardy-Weinberg Equilibrium exact test p-value threshold (default: 0)
- `Fchr` <STRING>: Filter specific chromosome (default: empty string for no filter)

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
