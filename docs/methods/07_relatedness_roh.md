# Relatedness Estimation and Runs of Homozygosity

**Input dataset**: Unlinked SNP dataset (`biallelic_snps_thinned.vcf.gz` and PLINK format)

---

## Relatedness Estimation

### Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[N_SNPS]` | SNP count in final analysis VCF |
| `[N_SAMPLES]` | Number of individuals |

### Methods text

Pairwise relatedness was estimated using four complementary estimators, each suited to different assumptions and reference populations:

**Yang et al. (2010) genomic relatedness** was computed using vcftools (Danecek et al. 2011) (`--relatedness`), which estimates the unadjusted A_jk statistic. Under this estimator, values near 1 indicate identity (monozygotic twins or technical duplicates) and values near 0 indicate unrelated individuals within the same population.

**Manichaikul et al. (2010) kinship** (the `--relatedness2` estimator in vcftools) was also computed; expected values are approximately 0.5 for monozygotic twins, 0.25 for first-degree relatives (parent–offspring or full siblings), 0.125 for second-degree relatives, and near 0 for unrelated individuals.

**PLINK IBD** (PI_HAT) was estimated using `plink --genome` (Chang et al. 2015) on the PLINK-format unlinked SNP dataset ([N_SNPS] SNPs, [N_SAMPLES] individuals), providing the proportion of the genome shared identical by descent.

**KING kinship** was calculated using `plink2 --make-king-table` on the same dataset, providing the robust KING-within estimator (Manichaikul et al. 2010) that is less sensitive to population structure than IBS-based estimators.

All four matrices were compared to assess consistency and to identify potential cryptic relatives or duplicate samples not removed during preprocessing.

---

## Runs of Homozygosity (ROH)

**Input dataset**: Unlinked SNP dataset (`biallelic_snps_thinned.vcf.gz`)

### Placeholders to fill in

| Placeholder | Where to find the value |
|---|---|
| `[ROH_GROUP_BY]` | `config["parameters"]["roh"]["group_by"]` (e.g. Site) |

### Methods text

Runs of homozygosity (ROH) were detected using the `roh` command in bcftools version [VERSION] (Danecek et al. 2021) on the unlinked SNP dataset. ROH length and number were summarised per individual and compared across [ROH_GROUP_BY] groups to assess inbreeding and effective population size differences. The fraction of the genome covered by ROH (F_ROH) was used as an individual-level inbreeding coefficient.

---

## Key software and citations

- **Yang relatedness** – Yang, J. et al. (2010). Common SNPs explain a large proportion of the heritability for human height. *Nature Genetics*, 42, 565–569. https://doi.org/10.1038/ng.608
- **Manichaikul / KING** – Manichaikul, A. et al. (2010). Robust relationship inference in genome-wide association studies. *Bioinformatics*, 26, 2867–2873. https://doi.org/10.1093/bioinformatics/btq559
- **vcftools** – Danecek, P. et al. (2011). The variant call format and VCFtools. *Bioinformatics*, 27, 2156–2158. https://doi.org/10.1093/bioinformatics/btr330
- **PLINK 2** – Chang, C.C. et al. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, 4, 7. https://doi.org/10.1186/s13742-015-0047-8
- **bcftools ROH** – Danecek, P. et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10, giab008. https://doi.org/10.1093/gigascience/giab008
