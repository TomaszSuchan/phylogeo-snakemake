# VCF Filtering and Dataset Preparation

## Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[N_SAMPLES_RAW]` | Total samples in the raw ipyrad VCF |
| `[N_SAMPLES_KEPT]` | Samples retained after all filtering |
| `[N_REMOVED_RELATED]` | Samples removed for relatedness (if filtering enabled) |
| `[KING_THRESHOLD]` | `config["parameters"]["relatedness_filtering"]["king_threshold"]` (default: 0.0884) |
| `[MAC_THRESHOLD]` | `config["parameters"]["relatedness_filtering"]["mac_threshold"]` (default: 1) |
| `[F_MISSING]` | `config["parameters"]["vcf_filtering"]["f_missing"]` (default: 1.0, i.e. disabled) |
| `[MAF]` | `config["parameters"]["vcf_filtering"]["maf"]` (default: 0, i.e. disabled) |
| `[N_SNPS_BIALLELIC]` | SNP count in `biallelic_snps.vcf.gz` (from VCF stats output) |
| `[N_RAD_LOCI]` | Number of RAD loci retained |
| `[N_SNPS_FINAL]` | SNP count in the final analysis VCF |
| `[THINNING_METHOD]` | `config["parameters"]["vcf_thinning"]["method"]` (default: `max_coverage`) |
| `[LD_R2]` | `config["parameters"]["ld_pruning"]["r2"]` (default: 0.5) |
| `[LD_WINDOW]` | `config["parameters"]["ld_pruning"]["window_size"]` (default: 50) |
| `[LD_STEP]` | `config["parameters"]["ld_pruning"]["step_size"]` (default: 5) |

---

## Methods text

### RAD-seq assembly and variant calling

Sequencing data were assembled and variants called using ipyrad version [VERSION] (Eaton & Overcast 2020) with [DESCRIBE ASSEMBLY SETTINGS]. This yielded an initial dataset of [N_SAMPLES_RAW] individuals and [N_LOCI_TOTAL] RAD loci.

### Sample filtering

The ipyrad VCF output was coordinate-sorted and bgzip-compressed using vcf-sort and bgzip (part of the VCFtools package; Danecek et al. 2011) and indexed with bcftools version [VERSION] (Danecek et al. 2021).

<!-- Include the following paragraph only if relatedness_filtering.enabled = true -->
**[If relatedness filtering was applied]** Prior to population genetic analysis, pairwise kinship was estimated among all [N_SAMPLES_RAW] individuals using the KING-robust statistic as implemented in PLINK 2 version [VERSION] (Chang et al. 2015), computed from biallelic variants with a minimum minor allele count of [MAC_THRESHOLD]. Individuals were greedily excluded until no pair exceeded a KING coefficient of [KING_THRESHOLD] (approximately second-degree relatedness). This step removed [N_REMOVED_RELATED] individual(s), leaving [N_SAMPLES_KEPT] individuals for downstream analysis.

### SNP filtering and construction of analysis datasets

Starting from the [N_SAMPLES_KEPT]-sample VCF, the following variant-level filters were applied sequentially using bcftools:

<!-- Include the F_MISSING line only if f_missing < 1.0 -->
1. Variants with a per-site missing-data rate exceeding [F_MISSING] were removed.
2. Only biallelic SNPs with a minor allele count greater than one (MAC > 1) were retained<!-- , with a minimum minor allele frequency (MAF) of [MAF] --> .

This produced the **biallelic SNP dataset** comprising [N_SNPS_BIALLELIC] SNPs distributed across [N_RAD_LOCI] RAD loci.

<!-- Choose ONE of the three paragraphs below based on thinning_strategy, delete the others -->

**[For thinning_strategy: "thinning" — Unlinked SNP dataset]**
Because SNPs within the same RAD locus are physically linked and thus non-independent, one SNP was retained per RAD locus using a custom Python script. The selection criterion was maximum sample coverage (i.e. the SNP with the fewest missing genotypes was preferred; ties were broken at random). This produced the **unlinked SNP dataset** of [N_SNPS_FINAL] SNPs from [N_RAD_LOCI] RAD loci (`biallelic_snps_thinned.vcf.gz`), which was used for all analyses that assume linkage equilibrium.

**[For thinning_strategy: "ld_pruning" — LD-pruned SNP dataset]**
SNPs were pruned for linkage disequilibrium (LD) using PLINK 1.9 (Chang et al. 2015) with a sliding-window approach (window size [LD_WINDOW] SNPs, step size [LD_STEP] SNPs, r² threshold [LD_R2]). This yielded the **LD-pruned SNP dataset** of [N_SNPS_FINAL] SNPs (`biallelic_snps_ld_pruned.vcf.gz`), used for all analyses assuming SNP independence.

**[For thinning_strategy: "none" — All-biallelic SNP dataset]**
No LD-based thinning was applied; the full set of [N_SNPS_BIALLELIC] biallelic SNPs was retained as the **all-biallelic SNP dataset** (`biallelic_snps_all.vcf.gz`).

### All-sites dataset (for pixy only)

To enable unbiased estimation of nucleotide diversity (π), absolute pairwise divergence (dXY), and FST using pixy (Korunes & Samuk 2021), invariant sites are required. An **all-sites dataset** was reconstructed from the ipyrad `.loci` file using a custom Python script that merged invariant positions with polymorphic sites in locus order, retaining the [N_SAMPLES_KEPT] individuals. The resulting VCF was sorted and bgzip-compressed (`allsites.vcf.gz`).

---

## Summary of output datasets

| Dataset | File | Use |
|---|---|---|
| Unlinked SNP dataset | `biallelic_snps_thinned.vcf.gz` (+ `.bed/.bim/.fam`, `.str`) | Structure, ADMIXTURE, PCA, genetic distances, relatedness |
| LD-pruned SNP dataset | `biallelic_snps_ld_pruned.vcf.gz` | Same as above, alternative to thinning |
| All-biallelic SNP dataset | `biallelic_snps_all.vcf.gz` | Same as above, without LD filtering |
| All-sites dataset | `allsites.vcf.gz` | pixy (π, FST, dXY) only |

---

## Key software and citations

- **ipyrad** – Eaton, D.A.R. & Overcast, I. (2020). ipyrad: Interactive assembly and analysis of RADseq datasets. *Bioinformatics*, 36, 2592–2594. https://doi.org/10.1093/bioinformatics/btrzz287
- **bcftools / VCFtools** – Danecek, P. et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10, giab008. https://doi.org/10.1093/gigascience/giab008; Danecek, P. et al. (2011). The variant call format and VCFtools. *Bioinformatics*, 27, 2156–2158. https://doi.org/10.1093/bioinformatics/btr330
- **PLINK 2** – Chang, C.C. et al. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, 4, 7. https://doi.org/10.1186/s13742-015-0047-8
- **KING** – Manichaikul, A. et al. (2010). Robust relationship inference in genome-wide association studies. *Bioinformatics*, 26, 2867–2873. https://doi.org/10.1093/bioinformatics/btq559
