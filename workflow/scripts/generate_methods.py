"""
Generate a manuscript-ready methods.md for one project.

Sources used (all filled in automatically):
  - snakemake.params.analyses   – which tools are enabled
  - snakemake.params.parameters – every config parameter
  - data/mainparams             – STRUCTURE BURNIN / NUMREPS
  - vcf_stats.txt               – variant count, RAD-loci count, sample count
  - workflow/envs/*.yaml        – software versions

Best-K values are NOT written to the methods (they are results).
"""

import glob
import os
import re
import textwrap


# ── helpers ──────────────────────────────────────────────────────────────────

def parse_mainparams(path):
    """Return {KEY: value} from a STRUCTURE mainparams file."""
    out = {}
    try:
        with open(path) as fh:
            for line in fh:
                m = re.match(r'\s*#define\s+(\w+)\s+(\S+)', line)
                if m:
                    out[m.group(1)] = m.group(2)
    except FileNotFoundError:
        pass
    return out


def parse_vcf_stats(path):
    """Return dict with 'variants', 'rad_fragments', 'samples' from vcf_stats.txt."""
    result = {}
    try:
        with open(path) as fh:
            for line in fh:
                for key, label in [("variants",      "Number of variants:"),
                                   ("rad_fragments",  "Number of RAD fragments:"),
                                   ("samples",        "Number of samples:")]:
                    if line.startswith(label):
                        result[key] = line.split(":", 1)[1].strip()
    except FileNotFoundError:
        pass
    return result


def parse_versions(envs_dir):
    """
    Scan all *.yaml files in envs_dir and return {package_name_lower: version}.
    Handles dependency lines like:
        - bcftools=1.22
        - plink=1.90b6.12
        - r-pophelper=2.3.1
        - python>=3.9          (skipped – no pinned version)
    """
    versions = {}
    for yf in sorted(glob.glob(os.path.join(envs_dir, "*.yaml"))):
        try:
            with open(yf) as fh:
                in_deps = False
                for line in fh:
                    stripped = line.strip()
                    if stripped == "dependencies:":
                        in_deps = True
                        continue
                    if in_deps:
                        if stripped.startswith("-"):
                            dep = stripped.lstrip("- ").strip()
                            # Only accept pinned exact versions (=), not ranges (>=, <=, ~=)
                            m = re.match(r'^([\w\-\.]+)=([^=><~\s\*]+)', dep)
                            if m:
                                name = m.group(1).lower()
                                ver  = m.group(2)
                                # Don't overwrite a version already found in an earlier file
                                if name not in versions:
                                    versions[name] = ver
        except Exception:
            pass
    return versions


def v(versions, key, fallback="[VERSION]"):
    """Look up a software version, return fallback string if not found."""
    return versions.get(key.lower(), fallback)


# ── gather inputs ─────────────────────────────────────────────────────────────

analyses = snakemake.params.analyses        # {tool: bool}
p        = snakemake.params.parameters      # project parameters dict
project  = snakemake.params.project
envs_dir = snakemake.params.envs_dir

mainparams = parse_mainparams(snakemake.input.mainparams)
vcf_stats  = parse_vcf_stats(snakemake.input.vcf_stats)
versions   = parse_versions(envs_dir)

n_snps     = vcf_stats.get("variants",      "[N_SNPS]")
n_loci     = vcf_stats.get("rad_fragments", "[N_RAD_LOCI]")
n_samples  = vcf_stats.get("samples",       "[N_SAMPLES]")

# ── dataset label from thinning strategy ─────────────────────────────────────

thinning_strat = p.get("thinning_strategy", "thinning")
DATASET = {
    "thinning":   ("unlinked biallelic SNP dataset",       "biallelic_snps_thinned.vcf.gz"),
    "ld_pruning": ("LD-pruned biallelic SNP dataset",      "biallelic_snps_ld_pruned.vcf.gz"),
    "none":       ("all-biallelic SNP dataset",            "biallelic_snps_all.vcf.gz"),
}
dataset_label, dataset_file = DATASET.get(thinning_strat, DATASET["thinning"])

# ── section builders ─────────────────────────────────────────────────────────

sections = []


# 1. Data filtering ────────────────────────────────────────────────────────────

vcf_filt  = p.get("vcf_filtering", {})
f_missing = vcf_filt.get("f_missing", 1.0)
maf       = vcf_filt.get("maf",       0.0)
rel_filt  = p.get("relatedness_filtering", {})
rel_on    = rel_filt.get("enabled",         False)
king_thr  = rel_filt.get("king_threshold",  0.0884)
mac_thr   = rel_filt.get("mac_threshold",   1)
vcf_thin  = p.get("vcf_thinning", {})
thin_meth = vcf_thin.get("method", "max_coverage")
ld        = p.get("ld_pruning", {})
ld_r2     = ld.get("r2",          0.5)
ld_win    = ld.get("window_size", 50)
ld_step   = ld.get("step_size",   5)

filt_parts = []

if rel_on:
    filt_parts.append(
        f"Prior to SNP filtering, pairwise kinship was estimated using the "
        f"KING-robust statistic in PLINK 2 version {v(versions,'plink2')} "
        f"(Chang et al. 2015) from biallelic variants with MAC ≥ {mac_thr}. "
        f"Samples were greedily excluded until no pair exceeded a KING coefficient "
        f"of {king_thr} (approximately second-degree relatedness)."
    )

miss_clause = (f"Variants with a per-site missing-data rate above {f_missing} were "
               f"excluded. " if f_missing < 1.0 else "")
maf_clause  = (f" and a minimum minor allele frequency (MAF) of {maf}"
               if maf > 0 else "")
filt_parts.append(
    f"{miss_clause}Variant calls were restricted to biallelic SNPs with a minor "
    f"allele count greater than one (MAC > 1){maf_clause} using bcftools version "
    f"{v(versions,'bcftools')} (Danecek et al. 2021)."
)

if thinning_strat == "thinning":
    thin_desc = {
        "max_coverage": "the SNP with the greatest number of called genotypes "
                        "(maximum sample coverage; ties broken at random)",
        "random":       "a randomly selected SNP",
        "weighted":     "a coverage-weighted randomly selected SNP",
    }.get(thin_meth, thin_meth)
    filt_parts.append(
        f"Because SNPs within the same RAD locus share identical flanking sequence "
        f"and are therefore physically linked, one SNP was retained per RAD locus "
        f"by selecting {thin_desc}. This produced the **unlinked biallelic SNP dataset** "
        f"({n_snps} SNPs across {n_loci} RAD loci; {n_samples} individuals), "
        f"used for all analyses that assume linkage equilibrium."
    )
elif thinning_strat == "ld_pruning":
    filt_parts.append(
        f"SNPs were pruned for linkage disequilibrium (LD) using PLINK version "
        f"{v(versions,'plink')} (Chang et al. 2015) with a sliding window of "
        f"{ld_win} SNPs, step size {ld_step} SNPs, and an r² threshold of "
        f"{ld_r2}. This produced the **LD-pruned biallelic SNP dataset** ({n_snps} SNPs, "
        f"{n_samples} individuals)."
    )
else:
    filt_parts.append(
        f"No LD-based thinning was applied; all {n_snps} biallelic SNPs were "
        f"retained as the **all-biallelic SNP dataset** ({n_samples} individuals)."
    )

if analyses.get("pixy", False):
    filt_parts.append(
        f"An **all-sites dataset** including both variable and invariant positions "
        f"was reconstructed from the ipyrad ˬloci file using a custom Python "
        f"script, required for unbiased estimation of nucleotide diversity and "
        f"divergence with pixy (Korunes & Samuk 2021)."
    )

sections.append(("Data filtering and dataset construction",
                 "\n\n".join(filt_parts)))


# 2. Population structure ──────────────────────────────────────────────────────

k_vals     = p.get("k_values", list(range(1, 10)))
k_min, k_max = min(k_vals), max(k_vals)
struct_parts = []

if analyses.get("structure", False):
    burnin  = mainparams.get("BURNIN",  "100000")
    numreps = mainparams.get("NUMREPS", "200000")
    reps    = p.get("structure", {}).get("replicates", 10)
    struct_parts.append(
        f"Bayesian clustering was performed with STRUCTURE version "
        f"{v(versions,'structure')} (Pritchard et al. 2000) on the {dataset_label} "
        f"({n_snps} SNPs, {n_samples} individuals) using the admixture model with "
        f"correlated allele frequencies (FREQSCORR = 1, NOADMIX = 0, "
        f"INFERALPHA = 1), a burn-in of {int(burnin):,} MCMC iterations, and "
        f"{int(numreps):,} sampling iterations. K = {k_min}–{k_max} was "
        f"tested with {reps} independent replicates each. Replicate runs were "
        f"aligned using the FullSearch algorithm in pophelper version "
        f"{v(versions,'r-pophelper')} (Francis 2017), and the optimal K was "
        f"selected by the ΔK method (Evanno et al. 2005)."
    )

if analyses.get("faststructure", False):
    tol   = p.get("faststructure", {}).get("tol",   "10e-6")
    prior = p.get("faststructure", {}).get("prior", "simple")
    struct_parts.append(
        f"Ancestry inference was additionally performed with fastStructure "
        f"version {v(versions,'faststructure')} (Raj et al. 2014) on the "
        f"{dataset_label} ({n_snps} SNPs, {n_samples} individuals) using the "
        f"{prior} prior and a convergence tolerance of {tol}. K = "
        f"{k_min}–{k_max} was each run once; the optimal K was identified "
        f"with the chooseK.py utility supplied with fastStructure."
    )

if analyses.get("admixture", False):
    struct_parts.append(
        f"Maximum-likelihood ancestry estimation was carried out with ADMIXTURE "
        f"version {v(versions,'admixture')} (Alexander et al. 2009) on the "
        f"{dataset_label} ({n_snps} SNPs, {n_samples} individuals). "
        f"K = {k_min}–{k_max} was tested; the optimal K was "
        f"selected as the value minimising five-fold cross-validation error."
    )

if analyses.get("dapc", False):
    dp    = p.get("dapc", {})
    n_pca = dp.get("n_pca",       50)
    n_da  = dp.get("n_da",        10)
    crit  = dp.get("criterion",   "diffNgroup")
    k_max_dapc = dp.get("max_n_clust", 10)
    struct_parts.append(
        f"Discriminant Analysis of Principal Components (DAPC; Jombart et al. 2010) "
        f"was performed in R using adegenet on the {dataset_label} "
        f"({n_snps} SNPs, {n_samples} individuals), retaining {n_pca} PCA axes "
        f"prior to discriminant analysis and {n_da} discriminant functions. "
        f"The optimal number of clusters (K = 1–{k_max_dapc}) was "
        f"determined by K-means cross-validation (30 replicates) using the "
        f"'{crit}' criterion."
    )

if analyses.get("construct", False):
    cp       = p.get("construct", {})
    n_iter   = cp.get("n_iterations", 10000)
    n_chains = cp.get("n_chains",     1)
    struct_parts.append(
        f"Spatially explicit population structure was modelled with conStruct "
        f"(Bradburd et al. 2018) on the {dataset_label} "
        f"({n_snps} SNPs, {n_samples} individuals), testing "
        f"K = {k_min}–{k_max} spatial layers with {n_chains} "
        f"MCMC chain(s) of {n_iter:,} iterations each."
    )

if struct_parts:
    sections.append(("Population structure and ancestry",
                     "\n\n".join(struct_parts)))


# 3. PCA / ordination ─────────────────────────────────────────────────────────

pca_parts = []

if analyses.get("pcaone", False):
    pp    = p.get("PCAone", {})
    n_pc  = pp.get("PCnum",      10)
    svd_m = pp.get("SVD_method",  3)
    SVD_NAMES = {0: "IRAM", 1: "single-pass randomised SVD",
                 2: "window-based randomised SVD", 3: "full SVD"}
    svd_desc = SVD_NAMES.get(int(svd_m), f"method {svd_m}")
    pca_parts.append(
        f"Principal component analysis (PCA) was performed with PCAone version "
        f"{v(versions,'pcaone')} (Zhang & Meisner 2023) on the {dataset_label} "
        f"({n_snps} SNPs, {n_samples} individuals), computing the {n_pc} leading "
        f"principal components via {svd_desc}."
    )

if analyses.get("pcaone_emu", False) or analyses.get("emupca", False):
    pp   = p.get("PCAone", {})
    n_pc = pp.get("PCnum", 10)
    pca_parts.append(
        f"PCA was additionally run in EMU mode (expectation–maximisation "
        f"imputation of missing genotypes) in PCAone version "
        f"{v(versions,'pcaone')}, retaining {n_pc} principal components."
    )

if analyses.get("pcoa", False):
    pca_parts.append(
        f"Principal coordinates analysis (PCoA) was performed on pairwise "
        f"genetic distance matrices (see below) using the cmdscale function in R."
    )

if pca_parts:
    sections.append(("Principal component analysis",
                     "\n\n".join(pca_parts)))


# 4. Genetic diversity ────────────────────────────────────────────────────────

div_parts = []

if analyses.get("pixy", False):
    pp     = p.get("pixy", {})
    stats  = pp.get("stats",                ["pi"])
    win    = pp.get("window_size",          10000)
    nboot  = pp.get("bootstrap_replicates", 1000)
    STAT_LABELS = {"pi": "nucleotide diversity (π)",
                   "fst": "FₛT",
                   "dxy": "dₛY"}
    stats_str = ", ".join(STAT_LABELS.get(s, s) for s in stats)
    verb = "was" if len(stats) == 1 else "were"
    div_parts.append(
        f"{stats_str.capitalize()} {verb} estimated using pixy version "
        f"{v(versions,'pixy')} (Korunes & Samuk 2021) on the all-sites dataset, "
        f"which includes both variable and invariant positions, in non-overlapping "
        f"windows of {win:,} bp. Population-level summaries and 95 % "
        f"confidence intervals were obtained by bootstrapping across windows "
        f"({nboot:,} replicates)."
    )

if analyses.get("amova", False):
    ap     = p.get("amova", {})
    strata = ap.get("strata", [])
    nperm  = ap.get("nperm",  999)
    strata_str = " and ".join(strata) if strata else "[hierarchical groupings]"
    div_parts.append(
        f"Analysis of molecular variance (AMOVA; Excoffier et al. 1992) was "
        f"performed in R using the poppr package (Kamvar et al. 2014) on the "
        f"{dataset_label}, partitioning genetic variance hierarchically among "
        f"{strata_str}. Significance was assessed with {nperm} random permutations."
    )

if div_parts:
    sections.append(("Genetic diversity and differentiation",
                     "\n\n".join(div_parts)))


# 5. Genetic distances and networks ───────────────────────────────────────────

dist_parts = []

if analyses.get("gen_dist", False):
    dist_parts.append(
        f"Four pairwise individual-level genetic distance matrices were computed "
        f"from the {dataset_label} ({n_snps} SNPs, {n_samples} individuals) using "
        f"custom Python scripts operating on PLINK binary files via the bed-reader "
        f"library: (i) Kosman–Leonard distance (Kosman & Leonard 2005), "
        f"a standardised distance for diploid codominant markers; "
        f"(ii) Euclidean distance from per-SNP dosage vectors (0, 1, 2) with "
        f"mean imputation of missing values; "
        f"(iii) p-distance (normalised Hamming distance), averaging "
        f"|gᵢ − gⱼ| / 2 over loci with data in "
        f"both individuals; and (iv) average squared difference "
        f"(bed2diffs-style Gram-matrix formulation)."
    )

if analyses.get("neighbornet", False):
    nn_color  = p.get("neighbornet", {}).get("color_by", ["Site"])
    color_str = ", ".join(nn_color) if isinstance(nn_color, list) else nn_color
    dist_parts.append(
        f"A NeighborNet split network (Bryant & Moulton 2004) was constructed "
        f"from the p-distance matrix using the phangorn package (Schliep 2011) "
        f"in R and visualised with tanggle, with tips coloured by {color_str}."
    )

if dist_parts:
    sections.append(("Genetic distances and networks",
                     "\n\n".join(dist_parts)))


# 6. Phylogenetics ────────────────────────────────────────────────────────────

phy_parts = []

for iq_key, label in [("iqtree",        ""),
                       ("iqtree_trimal", " with trimAl alignment preprocessing"),
                       ("iqtree_robust", " (robust-phylogeny mode)")]:
    if not analyses.get(iq_key, False):
        continue
    ip       = p.get("iqtree", {})
    model    = ip.get("model",             "MFP")
    boots    = ip.get("bootstraps",        1000)
    outgroup = ip.get("outgroup",          None)
    supp_thr = ip.get("support_threshold", 70)

    model_clause = (
        "the best-fit substitution model identified by ModelFinder "
        "(Kalyaanamoorthy et al. 2017) across the full model space"
        if model == "MFP" else f"the {model} substitution model"
    )
    root_clause = (f"rooted on {outgroup}" if outgroup else "midpoint-rooted")

    trimal_clause = ""
    if iq_key == "iqtree_trimal":
        gt = p.get("trimal", {}).get("gt", 0.95)
        trimal_clause = (
            f" Prior to tree inference, the alignment was trimmed with trimAl "
            f"version {v(versions,'trimal')} (Capella-Gutiérrez et al. 2009) "
            f"using a gap threshold of {gt} (columns with >{100*(1-gt):.0f} % "
            f"gaps removed)."
        )

    phy_parts.append(
        f"A maximum-likelihood phylogenetic tree{label} was inferred from the "
        f"concatenated RAD-seq alignment produced by ipyrad using IQ-TREE 2 "
        f"version {v(versions,'iqtree')} (Minh et al. 2020) with "
        f"{model_clause}.{trimal_clause} Branch support was assessed with "
        f"{boots:,} ultrafast bootstrap replicates (Hoang et al. 2018); nodes "
        f"below {supp_thr} % support are not shown. The tree was "
        f"{root_clause} for display."
    )

if analyses.get("fineradstructure", False):
    fr  = p.get("fineradstructure", {})
    cl  = fr.get("cluster", {})
    mc  = cl.get("mcmc_iterations", 1000000)
    bi  = cl.get("burnin",          1000000)
    th  = cl.get("thinning",        1000)
    mct = fr.get("tree", {}).get("mcmc_iterations", 10000)
    frs_ver = v(versions, "fineradstructure")
    phy_parts.append(
        f"Co-ancestry-based population clustering was performed with "
        f"fineRADstructure version {frs_ver} (Malinsky et al. 2018). RADpainter "
        f"computed a pairwise haplotype co-ancestry matrix from the biallelic "
        f"SNP dataset ({n_samples} individuals), exploiting phase information "
        f"within individual RAD loci. fineSTRUCTURE (Lawson et al. 2012) then "
        f"clustered individuals via MCMC ({mc:,} sampling iterations, {bi:,} "
        f"burn-in, thinning every {th:,} iterations). A maximum-clade-credibility "
        f"population tree was built from {mct:,} tree-building iterations."
    )

if phy_parts:
    sections.append(("Phylogenetic analysis", "\n\n".join(phy_parts)))


# 7. Relatedness and ROH ──────────────────────────────────────────────────────

rel_parts = []

if analyses.get("relatedness", False):
    rel_parts.append(
        f"Pairwise relatedness among {n_samples} individuals was estimated from "
        f"the {dataset_label} using four complementary estimators: "
        f"(i) the Aⱼₖ genomic relatedness statistic (Yang et al. 2010) "
        f"via vcftools version {v(versions,'vcftools')} (Danecek et al. 2011) "
        f"(––relatedness); "
        f"(ii) the Manichaikul et al. (2010) kinship estimator "
        f"(––relatedness2); "
        f"(iii) PLINK IBD (PI⁠_HAT) from plink version "
        f"{v(versions,'plink')} ––genome; and "
        f"(iv) the KING-robust kinship coefficient from plink2 "
        f"––make-king-table (Manichaikul et al. 2010)."
    )

if analyses.get("roh", False):
    roh_group = p.get("roh", {}).get("group_by", ["Site"])
    group_str = " and ".join(roh_group) if isinstance(roh_group, list) else roh_group
    rel_parts.append(
        f"Runs of homozygosity (ROH) were detected with bcftools roh version "
        f"{v(versions,'bcftools')} (Danecek et al. 2021) on the {dataset_label}. "
        f"ROH length and abundance were summarised per individual and compared "
        f"across {group_str} groups; the fraction of the genome within ROH "
        f"(FᴯⲋH) was used as an individual inbreeding coefficient."
    )

if rel_parts:
    sections.append(("Relatedness and inbreeding",
                     "\n\n".join(rel_parts)))


# 8. Spatial genetics ─────────────────────────────────────────────────────────

if analyses.get("mapi", False):
    mp        = p.get("mapi", {})
    n_perm    = mp.get("n_permutations", 1000)
    halfwidth = mp.get("grid_halfwidth", 100000)
    alpha     = mp.get("alpha",           0.05)
    crs_proj  = mp.get("crs_projected",   "EPSG:8857")
    mapi_ver  = v(versions, "r-mapi")
    sections.append((
        "Spatial genetic structure",
        f"Spatial patterns of genetic differentiation were explored with MAPI "
        f"version {mapi_ver} (Rezende et al. 2021) using the Euclidean genetic "
        f"distance matrix ({n_samples} individuals). Sample coordinates were "
        f"projected to {crs_proj} for distance calculations over a hexagonal grid "
        f"with a cell half-width of {halfwidth:,} m. Significance of upper "
        f"and lower tails of the cell-value distribution was assessed by "
        f"{n_perm:,} label permutations (α = {alpha})."
    ))


# ── references (only for used tools) ─────────────────────────────────────────

refs = {}

# Always present
refs["bcftools"] = (
    "Danecek, P. et al. (2021). Twelve years of SAMtools and BCFtools. "
    "*GigaScience*, 10, giab008. https://doi.org/10.1093/gigascience/giab008"
)
refs["plink"] = (
    "Chang, C.C. et al. (2015). Second-generation PLINK: rising to the challenge "
    "of larger and richer datasets. *GigaScience*, 4, 7. "
    "https://doi.org/10.1186/s13742-015-0047-8"
)

if analyses.get("structure", False):
    refs["structure"] = (
        "Pritchard, J.K., Stephens, M. & Donnelly, P. (2000). Inference of "
        "population structure using multilocus genotype data. *Genetics*, 155, "
        "945–959."
    )
    refs["pophelper"] = (
        "Francis, R.M. (2017). pophelper: an R package and web app to analyse and "
        "visualise population structure. *Molecular Ecology Resources*, 17, 27–32."
    )
    refs["evanno"] = (
        "Evanno, G., Regnaut, S. & Goudet, J. (2005). Detecting the number of "
        "clusters using STRUCTURE: a simulation study. *Molecular Ecology*, 14, "
        "2611–2620."
    )

if analyses.get("faststructure", False):
    refs["faststructure"] = (
        "Raj, A., Stephens, M. & Pritchard, J.K. (2014). fastSTRUCTURE: variational "
        "inference of population structure in large SNP datasets. *Genetics*, 197, "
        "573–589."
    )

if analyses.get("admixture", False):
    refs["admixture"] = (
        "Alexander, D.H., Novembre, J. & Lange, K. (2009). Fast model-based "
        "estimation of ancestry in unrelated individuals. *Genome Research*, 19, "
        "1655–1664."
    )

if analyses.get("dapc", False):
    refs["dapc"] = (
        "Jombart, T., Devillard, S. & Balloux, F. (2010). Discriminant analysis of "
        "principal components: a new method for the analysis of genetically structured "
        "populations. *BMC Genetics*, 11, 94."
    )

if analyses.get("construct", False):
    refs["construct"] = (
        "Bradburd, G.S., Coop, G.M. & Ralph, P.L. (2018). Inferring continuous and "
        "discrete population genetic structure across space. *Genetics*, 210, 33–52."
    )

if analyses.get("pcaone", False) or analyses.get("pcaone_emu", False) \
        or analyses.get("emupca", False):
    refs["pcaone"] = (
        "Zhang, C. & Meisner, J. (2023). Fast and accurate out-of-core PCA framework "
        "for large biobank data. *Genome Research*, 33, 1599–1608."
    )

if analyses.get("pixy", False):
    refs["pixy"] = (
        "Korunes, K.L. & Samuk, K. (2021). pixy: Unbiased estimation of nucleotide "
        "diversity and divergence in the presence of missing data. "
        "*Molecular Ecology Resources*, 21, 1359–1368."
    )

if analyses.get("amova", False):
    refs["amova"] = (
        "Excoffier, L., Smouse, P.E. & Quattro, J.M. (1992). Analysis of molecular "
        "variance inferred from metric distances among DNA haplotypes. "
        "*Genetics*, 131, 479–491."
    )
    refs["poppr"] = (
        "Kamvar, Z.N., Tabima, J.F. & Grünwald, N.J. (2014). Poppr: an R "
        "package for genetic analysis of populations with clonal, partially clonal, "
        "and/or sexual reproduction. *PeerJ*, 2, e281."
    )

if analyses.get("gen_dist", False):
    refs["kosman"] = (
        "Kosman, E. & Leonard, K.J. (2005). Similarity coefficients for molecular "
        "markers in studies of genetic relationships between individuals. "
        "*Molecular Ecology*, 14, 415–424."
    )

if analyses.get("neighbornet", False):
    refs["neighbornet"] = (
        "Bryant, D. & Moulton, V. (2004). Neighbor-Net: an agglomerative method "
        "for the construction of phylogenetic networks. "
        "*Molecular Biology and Evolution*, 21, 255–265."
    )
    refs["phangorn"] = (
        "Schliep, K.P. (2011). phangorn: phylogenetic analysis in R. "
        "*Bioinformatics*, 27, 592–593."
    )

if any(analyses.get(k, False) for k in ("iqtree", "iqtree_trimal", "iqtree_robust")):
    refs["iqtree"] = (
        "Minh, B.Q. et al. (2020). IQ-TREE 2: new models and methods for "
        "phylogenetic inference. *Molecular Biology and Evolution*, 37, 1530–1534."
    )
    refs["modelfinder"] = (
        "Kalyaanamoorthy, S. et al. (2017). ModelFinder: fast model selection for "
        "accurate phylogenetic estimates. *Nature Methods*, 14, 587–589."
    )
    refs["ufboot"] = (
        "Hoang, D.T. et al. (2018). UFBoot2: improving the ultrafast bootstrap "
        "approximation. *Molecular Biology and Evolution*, 35, 518–522."
    )

if analyses.get("iqtree_trimal", False):
    refs["trimal"] = (
        "Capella-Gutiérrez, S., Silla-Martínez, J.M. & Gabaldon, T. (2009). "
        "trimAl: a tool for automated alignment trimming in large-scale phylogenetic "
        "studies. *Bioinformatics*, 25, 1972–1973."
    )

if analyses.get("fineradstructure", False):
    refs["fineradstructure"] = (
        "Malinsky, M. et al. (2018). RADpainter and fineRADstructure: population "
        "inference from RADseq data. *Molecular Biology and Evolution*, 35, "
        "1284–1290."
    )
    refs["finestructure"] = (
        "Lawson, D.J. et al. (2012). Inference of population structure using dense "
        "haplotype data. *PLOS Genetics*, 8, e1002453."
    )

if analyses.get("relatedness", False):
    refs["yang2010"] = (
        "Yang, J. et al. (2010). Common SNPs explain a large proportion of the "
        "heritability for human height. *Nature Genetics*, 42, 565–569."
    )
    refs["manichaikul"] = (
        "Manichaikul, A. et al. (2010). Robust relationship inference in genome-wide "
        "association studies. *Bioinformatics*, 26, 2867–2873."
    )
    refs["vcftools"] = (
        "Danecek, P. et al. (2011). The variant call format and VCFtools. "
        "*Bioinformatics*, 27, 2156–2158."
    )

if analyses.get("mapi", False):
    refs["mapi"] = (
        "Rezende, F. et al. (2021). MAPI: an R package for mapping averaged pairwise "
        "information. *Methods in Ecology and Evolution*, 12, 1305–1310."
    )


# ── assemble document ─────────────────────────────────────────────────────────

lines = [f"# Methods — {project}", ""]

for title, body in sections:
    lines.append(f"## {title}")
    lines.append("")
    for para in body.split("\n\n"):
        lines.append(textwrap.fill(para, width=92,
                                   break_long_words=False,
                                   break_on_hyphens=False))
        lines.append("")

lines.append("## References")
lines.append("")
for ref in refs.values():
    lines.append(f"- {ref}")
lines.append("")

with open(snakemake.output.methods, "w") as fh:
    fh.write("\n".join(lines))
