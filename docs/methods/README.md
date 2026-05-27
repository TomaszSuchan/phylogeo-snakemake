# Methods Documentation

This directory contains manuscript-ready methods descriptions for each analysis tool in the phylogeographic pipeline. Each file is designed to be copied and pasted directly into the Methods section of a manuscript, with `[PLACEHOLDER]` markers indicating values that must be filled in from your specific dataset or `config/config.yaml`.

## Dataset Naming Conventions

The pipeline produces several filtered datasets, referred to by the following human-readable names throughout all methods files. The corresponding file names under `results/<project>/filtered_data/` are:

| Human-readable name | File suffix | Description |
|---|---|---|
| **Unlinked SNP dataset** | `biallelic_snps_thinned.vcf.gz` / `.bed/.bim/.fam` / `.str` | Biallelic SNPs with one SNP selected per RAD locus (thinning strategy). Used for ancestry, PCA, genetic distance, and most population genetic analyses where SNP independence is required. |
| **LD-pruned SNP dataset** | `biallelic_snps_ld_pruned.vcf.gz` | Biallelic SNPs after PLINK-based LD pruning across all loci. Alternative to per-locus thinning. |
| **All-biallelic SNP dataset** | `biallelic_snps_all.vcf.gz` | All biallelic SNPs with MAC > 1, without LD filtering (thinning_strategy: none). |
| **All-sites dataset** | `allsites.vcf.gz` | All variable and invariant sites reconstructed from the ipyrad `.loci` file. Required for unbiased estimation of π, dXY, and FST with pixy. |

> The active dataset used by most analyses depends on `thinning_strategy` in `config.yaml`. With the default (`"thinning"`), analyses use the **unlinked SNP dataset**. The pipeline exports this to PLINK (`.bed/.bim/.fam`) and STRUCTURE (`.str`) formats regardless of which strategy is active.

## Available Methods Files

| File | Tools covered |
|---|---|
| [01_filtering.md](01_filtering.md) | VCF preprocessing, sample filtering, SNP selection, dataset generation |
| [02_population_structure.md](02_population_structure.md) | STRUCTURE, fastStructure, ADMIXTURE, DAPC, conStruct |
| [03_pca_ordination.md](03_pca_ordination.md) | PCAone, EMU-PCA, PCoA |
| [04_genetic_diversity.md](04_genetic_diversity.md) | pixy (π, FST, dXY), AMOVA |
| [05_genetic_distances_networks.md](05_genetic_distances_networks.md) | Kosman distance, Euclidean distance, p-distance, NeighborNet |
| [06_phylogenetics.md](06_phylogenetics.md) | IQ-TREE, trimAl + IQ-TREE, fineRADstructure |
| [07_relatedness_roh.md](07_relatedness_roh.md) | Yang relatedness, Manichaikul relatedness, PLINK IBD, KING, ROH |
| [08_spatial_genetics.md](08_spatial_genetics.md) | MAPI |

## How to use

1. Open the relevant `.md` file for the analysis you are describing.
2. Copy the methods paragraph(s) you need.
3. Replace all `[PLACEHOLDER]` tokens with values from your dataset or `config/config.yaml`.
4. Adjust wording as needed for your manuscript's style.
5. Add the cited references to your bibliography.
