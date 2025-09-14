# phylogeo-snakemake

This repository provides a **Snakemake pipeline** for population genetics and phylogeographic data analysis.  
It is designed to work directly with **ipyrad** output files — currently making use of the `.vcf` and `.loci` files (the latter is required to retrieve invariant sites for genetic diversity calculations).

All configuration is handled through `config/config.yaml`.

---

## Workflow

### Filtering
- Sort, compress, and index the VCF file (ipyrad output is often unsorted).  
- Filter for **biallelic SNPs**.  
- Thin the VCF to retain **one SNP per RAD fragment**.  
  - This step is customizable via `config.yaml`:  
    - Filter by coverage.  
    - Choose how SNPs are selected (randomly, by maximum coverage, or randomly weighted by coverage).  
    - If multiple SNPs have equal maximum coverage, the pipeline can select the first or a random SNP (default).  
- Export the thinned VCF to **PLINK** and **Structure** formats.  
- Apply additional filtering for varying missing data thresholds and export the resulting datasets to PLINK (for PCA analyses).  

### Analyses
- **PCA (PCAone):** run on datasets with different missing data thresholds, with missing data handled using the `--emu` method.  
- **PCA (VCF2PCACluster):** allows flexible SNP filtering, testing different MAF and missing data thresholds, and accounting for kinship using multiple methods.  
- **Structure analysis:** for a user-defined range of *K* and replicates. Runs for each *K* are plotted to assess consistency and convergence.  
- **Admixture analysis:** for a user-defined range of *K*.  
- **fastStructure analysis:** for a user-defined range of *K*.  
- **Relatedness analysis:** using two methods implemented in `vcftools`.  
- **Genetic diversity metrics (pixy):** computes π, FST, and dXY using both variable and invariant sites (ensuring unbiased estimates).  