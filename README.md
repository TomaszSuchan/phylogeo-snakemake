# phylogeo-snakemake

This repository contains a Snakemake pipeline for population genetics and phylogeographic data analysis.

## Pipeline Steps, Tools, and Parameters

### 1. Preprocessing and Filtering

**Rule file:** [`preprocessing.smk`](workflow/rules/preprocessing.smk)  
**Input:** ipyrad output - defined by prefix, in `config/config.yaml` eg.: `test_data/Belenois_aurota_rmoutgroup`  
**Key parameters:**  
- Only biallelic SNPs with minor allele count (MAC) > 1 are kept.
- VCF thinning: one SNP per RAD fragment, using method and parameters from `config/config.yaml`:
  - `min_coverage`: Minimum number of samples with data (default: 0)
  - `method`: SNP selection method per fragment (`max_coverage`, `random`, or `weighted`)
  - `ties`: How to resolve ties whan SNPs have teh same coverage (`random` or `first`)
  - `ns_tag`: INFO field tag for sample count (default: `NS`)
  - `id_pattern`: Regex for RAD fragment ID (default: `loc(\d+)_`)

### 2. Principal Component Analysis (PCA)

**Rule file:** [`pcaone.smk`](workflow/rules/pcaone.smk)  
**Tool:** [PCAone](https://github.com/GenomicsPL/pcaone)  
**Parameters:**  
- SVD method selected via `SVD_method` in [`config.yaml`](config/config.yaml):
  - `0`: IRAM
  - `1`: Yu's single-pass Randomized SVD
  - `2`: PCAone window-based Randomized SVD
  - `3`: Full SVD (default)
- Output:  
  - `pcaone_EMU/PCA_EMU.eigvecs`, `pcaone_EMU/PCA_EMU.eigvals`
  - `pcaone/PCA.eigvecs`, `pcaone/PCA.eigvals`

### 3. Population Structure Analysis – fastStructure

**Rule file:** [`faststructure.smk`](workflow/rules/faststructure.smk)  
**Tool:** [fastStructure](https://github.com/rajanil/fastStructure)  
**Parameters from config:**  
- `k_values`: List of K (number of clusters) to test (e.g., `[2, 3]`)
- `faststructure_tol`: Convergence criterion (default: `10e-6`)
- `faststructure_prior`: Prior choice (`simple` or `logistic`, default: `simple`)
- Output:  
  - `faststructure/input_K{k}.meanQ`, `faststructure/input_K{k}.meanP`
  - `faststructure/chooseK_results.txt`

### 4. Population Structure Analysis – ADMIXTURE

**Rule file:** [`admixture.smk`](workflow/rules/admixture.smk)  
**Tool:** [ADMIXTURE](https://www.genetics.ucla.edu/software/admixture/)  
**Parameters:**  
- `k_values`: List of K (number of clusters) to test (e.g., `[2, 3]`)
- Output:  
  - `admixture/input.{k}.Q`, `admixture/input.{k}.P`

## How to Run

1. Edit [`config/config.yaml`](config/config.yaml) to set your input and parameters.
2. Run:
   ```sh
   snakemake --use-conda --cores <number_of_cores>
   ```

To run on an HPC cluster with SLURM:

```sh
conda install bioconda::snakemake-executor-plugin-slurm
snakemake --executor slurm --default-resources --jobs N --use-conda
```

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- [Conda](https://docs.conda.io/)
- Tools: PCAone, fastStructure, ADMIXTURE, bcftools, vcftools (see rule files for details)

---