"""
Generate a manuscript-ready methods.md for one project.

Reads all analysis parameters from snakemake.params, the project config, and
the STRUCTURE mainparams file. Only includes sections for analyses that are
enabled in the project config. Actual sample count is read from the
samples_to_keep.txt produced by the preprocessing pipeline.

Output: results/{project}/methods.md
"""

import re
import textwrap


# ---------------------------------------------------------------------------
# Helper: parse STRUCTURE mainparams
# ---------------------------------------------------------------------------

def parse_mainparams(path):
    """Return dict of key→value from a STRUCTURE mainparams file."""
    values = {}
    try:
        with open(path) as fh:
            for line in fh:
                m = re.match(r'\s*#define\s+(\w+)\s+(\S+)', line)
                if m:
                    values[m.group(1)] = m.group(2)
    except FileNotFoundError:
        pass
    return values


# ---------------------------------------------------------------------------
# Helper: read sample count
# ---------------------------------------------------------------------------

def count_samples(path):
    try:
        with open(path) as fh:
            n = sum(1 for line in fh if line.strip())
        return str(n)
    except FileNotFoundError:
        return "[N_SAMPLES]"


# ---------------------------------------------------------------------------
# Gather inputs
# ---------------------------------------------------------------------------

analyses = snakemake.params.analyses          # dict: tool → bool
p        = snakemake.params.parameters        # project parameters dict
project  = snakemake.params.project

mainparams  = parse_mainparams(snakemake.input.mainparams)
n_samples   = count_samples(snakemake.input.samples_to_keep)

# Convenience accessors
k_vals          = p.get("k_values", list(range(1, 10)))
k_min, k_max    = min(k_vals), max(k_vals)
thinning_strat  = p.get("thinning_strategy", "thinning")
vcf_filt        = p.get("vcf_filtering", {})
f_missing       = vcf_filt.get("f_missing", 1.0)
maf             = vcf_filt.get("maf", 0.0)
rel_filt        = p.get("relatedness_filtering", {})
rel_enabled     = rel_filt.get("enabled", False)
king_thresh     = rel_filt.get("king_threshold", 0.0884)
mac_thresh      = rel_filt.get("mac_threshold", 1)
vcf_thin        = p.get("vcf_thinning", {})
thin_method     = vcf_thin.get("method", "max_coverage")
ld              = p.get("ld_pruning", {})
ld_r2           = ld.get("r2", 0.5)
ld_win          = ld.get("window_size", 50)
ld_step         = ld.get("step_size", 5)

# Dataset name based on strategy
DATASET_NAMES = {
    "thinning":   ("unlinked SNP dataset",       "biallelic_snps_thinned.vcf.gz"),
    "ld_pruning": ("LD-pruned SNP dataset",       "biallelic_snps_ld_pruned.vcf.gz"),
    "none":       ("all-biallelic SNP dataset",   "biallelic_snps_all.vcf.gz"),
}
dataset_label, dataset_file = DATASET_NAMES.get(thinning_strat, DATASET_NAMES["thinning"])

# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

sections = []

# ── 1. Data filtering ────────────────────────────────────────────────────────

filt_parts = []

# Relatedness
if rel_enabled:
    filt_parts.append(
        f"Prior to SNP filtering, pairwise kinship was estimated using the "
        f"KING-robust statistic in PLINK 2 (Chang et al. 2015) from biallelic "
        f"variants with MAC ≥ {mac_thresh}. Samples were greedily excluded until "
        f"no pair exceeded a KING coefficient of {king_thresh} (approximately "
        f"second-degree relatedness)."
    )

# Missing-data filter
if f_missing < 1.0:
    filt_parts.append(
        f"Variants with a per-site missing-data rate above {f_missing} were removed."
    )

# Biallelic / MAF filter
maf_clause = f" and a minimum minor allele frequency (MAF) of {maf}" if maf > 0 else ""
filt_parts.append(
    f"Variant calls were restricted to biallelic SNPs with a minor allele "
    f"count greater than one (MAC > 1){maf_clause}, using bcftools (Danecek et al. 2021)."
)

# Thinning / LD
if thinning_strat == "thinning":
    method_desc = {
        "max_coverage": "the SNP with the fewest missing genotypes (maximum sample coverage)",
        "random": "a randomly selected SNP",
        "weighted": "a coverage-weighted randomly selected SNP",
    }.get(thin_method, thin_method)
    filt_parts.append(
        f"Because SNPs within the same RAD locus are physically linked, one SNP "
        f"per RAD locus was retained by selecting {method_desc} (ties broken at "
        f"random). This produced the **unlinked SNP dataset** "
        f"(`{dataset_file}`), used for all analyses that assume linkage equilibrium."
    )
elif thinning_strat == "ld_pruning":
    filt_parts.append(
        f"SNPs were pruned for linkage disequilibrium using PLINK 1.9 "
        f"(Chang et al. 2015) with a sliding window of {ld_win} SNPs, step size "
        f"{ld_step} SNPs, and an r² threshold of {ld_r2}. This produced the "
        f"**LD-pruned SNP dataset** (`{dataset_file}`)."
    )
else:
    filt_parts.append(
        f"No LD-based thinning was applied; all biallelic SNPs were retained "
        f"as the **all-biallelic SNP dataset** (`{dataset_file}`)."
    )

# All-sites for pixy
if analyses.get("pixy", False):
    filt_parts.append(
        "An **all-sites dataset** including both variant and invariant positions "
        "was reconstructed from the ipyrad `.loci` file using a custom Python "
        "script, required for unbiased estimation of nucleotide diversity and "
        "divergence with pixy (Korunes & Samuk 2021)."
    )

sections.append(("Data filtering and dataset construction",
                 "\n\n".join(filt_parts)))

# ── 2. Population structure ──────────────────────────────────────────────────

struct_parts = []

if analyses.get("structure", False):
    burnin  = mainparams.get("BURNIN",  "100000")
    numreps = mainparams.get("NUMREPS", "200000")
    reps    = p.get("structure", {}).get("replicates", 10)
    struct_parts.append(
        f"Bayesian clustering was performed with STRUCTURE version 2.3.4 "
        f"(Pritchard et al. 2000) on the {dataset_label} "
        f"({n_samples} individuals) using the admixture model with correlated "
        f"allele frequencies (FREQSCORR = 1, NOADMIX = 0, INFERALPHA = 1), "
        f"a burn-in of {int(burnin):,} MCMC iterations, and {int(numreps):,} "
        f"sampling iterations. K = {k_min}–{k_max} was tested with {reps} "
        f"independent replicates each. Replicate runs were aligned with the "
        f"FullSearch algorithm in pophelper v2.3.1 (Francis 2017), and the "
        f"optimal K was selected by the ΔK method (Evanno et al. 2005)."
    )

if analyses.get("faststructure", False):
    tol   = p.get("faststructure", {}).get("tol",   "10e-6")
    prior = p.get("faststructure", {}).get("prior", "simple")
    struct_parts.append(
        f"Ancestry inference was also performed with fastStructure version 1.0 "
        f"(Raj et al. 2014) on the {dataset_label} "
        f"({n_samples} individuals) using the {prior} prior and a convergence "
        f"tolerance of {tol}. K = {k_min}–{k_max} was each run once; the optimal "
        f"K was identified with the chooseK.py utility."
    )

if analyses.get("admixture", False):
    struct_parts.append(
        f"Maximum-likelihood ancestry estimation was carried out with ADMIXTURE "
        f"version 1.3.0 (Alexander et al. 2009) on the {dataset_label} "
        f"({n_samples} individuals). K = {k_min}–{k_max} was tested; the optimal "
        f"K was selected as the value minimising five-fold cross-validation error."
    )

if analyses.get("dapc", False):
    dp = p.get("dapc", {})
    n_pca  = dp.get("n_pca",  50)
    n_da   = dp.get("n_da",   10)
    crit   = dp.get("criterion", "diffNgroup")
    k_dapc = dp.get("max_n_clust", 10)
    struct_parts.append(
        f"Discriminant Analysis of Principal Components (DAPC; Jombart et al. 2010) "
        f"was performed in R using adegenet on the {dataset_label} "
        f"({n_samples} individuals), retaining {n_pca} PCA axes prior to "
        f"discriminant analysis and {n_da} discriminant functions. The optimal "
        f"number of clusters (K = 1–{k_dapc}) was determined by K-means "
        f"cross-validation (30 replicates) using the '{crit}' criterion."
    )

if analyses.get("construct", False):
    cp = p.get("construct", {})
    n_iter   = cp.get("n_iterations", 10000)
    n_chains = cp.get("n_chains",     1)
    struct_parts.append(
        f"Spatially explicit population structure was modelled with conStruct "
        f"version 1.0.4 (Bradburd et al. 2018) on the {dataset_label}, "
        f"testing K = {k_min}–{k_max} spatial layers with {n_chains} MCMC "
        f"chain(s) of {n_iter:,} iterations each."
    )

if struct_parts:
    sections.append(("Population structure and ancestry", "\n\n".join(struct_parts)))

# ── 3. PCA / ordination ──────────────────────────────────────────────────────

pca_parts = []

if analyses.get("pcaone", False):
    pp      = p.get("PCAone", {})
    n_pc    = pp.get("PCnum",       10)
    svd_m   = pp.get("SVD_method",   3)
    svd_names = {0: "IRAM", 1: "single-pass randomised SVD",
                 2: "window-based randomised SVD", 3: "full SVD"}
    svd_desc  = svd_names.get(int(svd_m), f"method {svd_m}")
    pca_parts.append(
        f"Principal component analysis was performed with PCAone (Zhang & Meisner 2023) "
        f"on the {dataset_label} ({n_samples} individuals) computing the "
        f"{n_pc} leading principal components via {svd_desc}."
    )

if analyses.get("pcaone_emu", False) or analyses.get("emupca", False):
    pp   = p.get("PCAone", {})
    n_pc = pp.get("PCnum", 10)
    pca_parts.append(
        f"PCA was additionally run in EMU mode (expectation–maximisation for "
        f"missing data) in PCAone, retaining {n_pc} principal components."
    )

if analyses.get("pcoa", False):
    pca_parts.append(
        f"Principal coordinates analysis (PCoA) was performed on pairwise "
        f"genetic distance matrices (see Genetic distances below) using "
        f"the cmdscale function in R."
    )

if pca_parts:
    sections.append(("Principal component analysis", "\n\n".join(pca_parts)))

# ── 4. Genetic diversity ─────────────────────────────────────────────────────

div_parts = []

if analyses.get("pixy", False):
    pp      = p.get("pixy", {})
    stats   = pp.get("stats",                ["pi"])
    win     = pp.get("window_size",          10000)
    nboot   = pp.get("bootstrap_replicates", 1000)
    stats_str = ", ".join(
        {"pi": "nucleotide diversity (π)", "fst": "FST", "dxy": "dXY"}.get(s, s)
        for s in stats
    )
    div_parts.append(
        f"{stats_str.capitalize()} {'was' if len(stats)==1 else 'were'} estimated "
        f"using pixy version [VERSION] (Korunes & Samuk 2021) on the all-sites dataset, "
        f"which includes both variant and invariant positions, in non-overlapping "
        f"windows of {win:,} bp. Population-level summaries and 95% confidence "
        f"intervals were obtained by bootstrapping across windows ({nboot:,} replicates)."
    )

if analyses.get("amova", False):
    ap = p.get("amova", {})
    strata = ap.get("strata", [])
    nperm  = ap.get("nperm",  999)
    strata_str = ", ".join(strata) if strata else "[strata]"
    div_parts.append(
        f"Analysis of molecular variance (AMOVA; Excoffier et al. 1992) was "
        f"performed in R using the poppr package (Kamvar et al. 2014) on the "
        f"{dataset_label}, partitioning genetic variance hierarchically among "
        f"{strata_str}. Significance was assessed with {nperm} random permutations."
    )

if div_parts:
    sections.append(("Genetic diversity and differentiation", "\n\n".join(div_parts)))

# ── 5. Genetic distances and networks ────────────────────────────────────────

dist_parts = []

if analyses.get("gen_dist", False):
    dist_parts.append(
        f"Four pairwise genetic distance matrices were computed from the "
        f"{dataset_label} ({n_samples} individuals) using custom Python scripts "
        f"via the bed-reader library: "
        f"(i) Kosman–Leonard distance (Kosman & Leonard 2005), a standardised "
        f"distance for diploid codominant markers; "
        f"(ii) Euclidean distance from per-SNP dosage vectors (0, 1, 2) with "
        f"mean imputation of missing values; "
        f"(iii) p-distance (normalised Hamming distance), averaging |g_i − g_j| / 2 "
        f"over loci with data in both individuals; and "
        f"(iv) average squared difference (bed2diffs-style Gram-matrix formulation)."
    )

if analyses.get("neighbornet", False):
    nn_color = p.get("neighbornet", {}).get("color_by", ["Site"])
    color_str = ", ".join(nn_color) if isinstance(nn_color, list) else nn_color
    dist_parts.append(
        f"A NeighborNet split network (Bryant & Moulton 2004) was constructed "
        f"from the p-distance matrix using phangorn (Schliep 2011) and visualised "
        f"with tanggle in R, with tips coloured by {color_str}."
    )

if dist_parts:
    sections.append(("Genetic distances and networks", "\n\n".join(dist_parts)))

# ── 6. Phylogenetics ─────────────────────────────────────────────────────────

phy_parts = []

for iq_key, label in [("iqtree", ""), ("iqtree_trimal", " with trimAl preprocessing"),
                      ("iqtree_robust", " (robust mode)")]:
    if not analyses.get(iq_key, False):
        continue
    ip       = p.get("iqtree", {})
    model    = ip.get("model", "MFP")
    boots    = ip.get("bootstraps", 1000)
    outgroup = ip.get("outgroup", None)
    supp_thr = ip.get("support_threshold", 70)

    model_clause = (
        "the best-fit substitution model identified by ModelFinder (Kalyaanamoorthy "
        "et al. 2017) across the full model space"
        if model == "MFP" else f"the {model} substitution model"
    )
    root_clause = (
        f"rooted on {outgroup}" if outgroup else "midpoint-rooted"
    )
    trimal_clause = ""
    if iq_key == "iqtree_trimal":
        gt = p.get("trimal", {}).get("gt", 0.95)
        trimal_clause = (
            f" Prior to tree inference, the alignment was trimmed with trimAl "
            f"(Capella-Gutiérrez et al. 2009) using a gap threshold of {gt} "
            f"(columns with >{100*(1-gt):.0f}% gaps removed)."
        )

    phy_parts.append(
        f"A maximum-likelihood phylogenetic tree{label} was inferred from the "
        f"concatenated RAD-seq alignment produced by ipyrad using IQ-TREE 2 "
        f"(Minh et al. 2020) with {model_clause}.{trimal_clause} Branch support "
        f"was assessed with {boots:,} ultrafast bootstrap replicates (Hoang et al. 2018); "
        f"nodes below {supp_thr}% support are not shown. The tree was {root_clause} "
        f"for display."
    )

if analyses.get("fineradstructure", False):
    fr = p.get("fineradstructure", {})
    cl = fr.get("cluster", {})
    mcmc_clust = cl.get("mcmc_iterations", 1000000)
    burnin_fr  = cl.get("burnin",          1000000)
    thin_fr    = cl.get("thinning",        1000)
    mcmc_tree  = fr.get("tree", {}).get("mcmc_iterations", 10000)
    phy_parts.append(
        f"Co-ancestry-based clustering was performed with fineRADstructure "
        f"(Malinsky et al. 2018). RADpainter computed pairwise haplotype "
        f"co-ancestry from the biallelic SNP dataset ({n_samples} individuals), "
        f"and fineSTRUCTURE (Lawson et al. 2012) clustered individuals via MCMC "
        f"({mcmc_clust:,} sampling iterations, {burnin_fr:,} burn-in, thinning "
        f"every {thin_fr:,} iterations). A maximum-clade-credibility tree was "
        f"built from {mcmc_tree:,} tree-building iterations."
    )

if phy_parts:
    sections.append(("Phylogenetic analysis", "\n\n".join(phy_parts)))

# ── 7. Relatedness / ROH ─────────────────────────────────────────────────────

rel_parts = []

if analyses.get("relatedness", False):
    rel_parts.append(
        f"Pairwise relatedness among {n_samples} individuals was estimated from "
        f"the {dataset_label} using four estimators: "
        f"(i) the A_jk genomic relatedness statistic (Yang et al. 2010) via "
        f"vcftools `--relatedness`; "
        f"(ii) the Manichaikul et al. (2010) kinship estimator via vcftools "
        f"`--relatedness2`; "
        f"(iii) PLINK IBD (PI_HAT) from `plink --genome`; and "
        f"(iv) the KING-robust kinship coefficient from `plink2 --make-king-table` "
        f"(Manichaikul et al. 2010)."
    )

if analyses.get("roh", False):
    roh_group = p.get("roh", {}).get("group_by", ["Site"])
    group_str = ", ".join(roh_group) if isinstance(roh_group, list) else roh_group
    rel_parts.append(
        f"Runs of homozygosity (ROH) were detected with `bcftools roh` "
        f"(Danecek et al. 2021) on the {dataset_label}. ROH length and abundance "
        f"were summarised per individual and compared across {group_str} groups; "
        f"the fraction of the genome in ROH (F_ROH) was used as an individual "
        f"inbreeding coefficient."
    )

if rel_parts:
    sections.append(("Relatedness and inbreeding", "\n\n".join(rel_parts)))

# ── 8. Spatial genetics ──────────────────────────────────────────────────────

if analyses.get("mapi", False):
    mp       = p.get("mapi", {})
    n_perm   = mp.get("n_permutations", 1000)
    halfwidth = mp.get("grid_halfwidth", 100000)
    alpha    = mp.get("alpha",           0.05)
    crs_proj = mp.get("crs_projected",   "EPSG:8857")
    sections.append((
        "Spatial genetic structure",
        f"Spatial patterns of genetic differentiation were explored with MAPI "
        f"(Rezende et al. 2021) using the Euclidean genetic distance matrix "
        f"({n_samples} individuals). Sample coordinates were projected to "
        f"{crs_proj} for distance calculations over a hexagonal grid with a "
        f"cell half-width of {halfwidth:,} m. Significance of upper and lower "
        f"tails of the cell-value distribution was assessed by {n_perm:,} "
        f"label-permutations (α = {alpha})."
    ))

# ---------------------------------------------------------------------------
# References (only for used tools)
# ---------------------------------------------------------------------------

refs = {}

refs["bcftools"] = (
    "Danecek, P. et al. (2021). Twelve years of SAMtools and BCFtools. "
    "*GigaScience*, 10, giab008. https://doi.org/10.1093/gigascience/giab008"
)
refs["vcftools"] = (
    "Danecek, P. et al. (2011). The variant call format and VCFtools. "
    "*Bioinformatics*, 27, 2156–2158. https://doi.org/10.1093/bioinformatics/btr330"
)
refs["plink"] = (
    "Chang, C.C. et al. (2015). Second-generation PLINK. "
    "*GigaScience*, 4, 7. https://doi.org/10.1186/s13742-015-0047-8"
)
if analyses.get("structure", False):
    refs["structure"] = (
        "Pritchard, J.K., Stephens, M. & Donnelly, P. (2000). Inference of "
        "population structure using multilocus genotype data. *Genetics*, 155, 945–959."
    )
    refs["pophelper"] = (
        "Francis, R.M. (2017). pophelper: an R package and web app to analyse and "
        "visualise population structure. *Molecular Ecology Resources*, 17, 27–32."
    )
    refs["evanno"] = (
        "Evanno, G., Regnaut, S. & Goudet, J. (2005). Detecting the number of "
        "clusters using STRUCTURE. *Molecular Ecology*, 14, 2611–2620."
    )
if analyses.get("faststructure", False):
    refs["faststructure"] = (
        "Raj, A., Stephens, M. & Pritchard, J.K. (2014). fastSTRUCTURE: variational "
        "inference of population structure in large SNP datasets. *Genetics*, 197, 573–589."
    )
if analyses.get("admixture", False):
    refs["admixture"] = (
        "Alexander, D.H., Novembre, J. & Lange, K. (2009). Fast model-based "
        "estimation of ancestry in unrelated individuals. "
        "*Genome Research*, 19, 1655–1664."
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
if analyses.get("pcaone", False) or analyses.get("pcaone_emu", False) or analyses.get("emupca", False):
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
        "variance inferred from metric distances among DNA haplotypes. *Genetics*, 131, 479–491."
    )
    refs["poppr"] = (
        "Kamvar, Z.N., Tabima, J.F. & Grünwald, N.J. (2014). Poppr: an R package for "
        "genetic analysis of populations. *PeerJ*, 2, e281."
    )
if analyses.get("gen_dist", False):
    refs["kosman"] = (
        "Kosman, E. & Leonard, K.J. (2005). Similarity coefficients for molecular "
        "markers in studies of genetic relationships. *Molecular Ecology*, 14, 415–424."
    )
if analyses.get("neighbornet", False):
    refs["neighbornet"] = (
        "Bryant, D. & Moulton, V. (2004). Neighbor-Net: an agglomerative method for "
        "the construction of phylogenetic networks. "
        "*Molecular Biology and Evolution*, 21, 255–265."
    )
    refs["phangorn"] = (
        "Schliep, K.P. (2011). phangorn: phylogenetic analysis in R. "
        "*Bioinformatics*, 27, 592–593."
    )
if any(analyses.get(k, False) for k in ("iqtree", "iqtree_trimal", "iqtree_robust")):
    refs["iqtree"] = (
        "Minh, B.Q. et al. (2020). IQ-TREE 2: new models and methods for phylogenetic "
        "inference. *Molecular Biology and Evolution*, 37, 1530–1534."
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
        "Capella-Gutiérrez, S., Silla-Martínez, J.M. & Gabaldón, T. (2009). trimAl: "
        "a tool for automated alignment trimming. *Bioinformatics*, 25, 1972–1973."
    )
if analyses.get("fineradstructure", False):
    refs["fineradstructure"] = (
        "Malinsky, M. et al. (2018). RADpainter and fineRADstructure: population "
        "inference from RADseq data. *Molecular Biology and Evolution*, 35, 1284–1290."
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
if analyses.get("mapi", False):
    refs["mapi"] = (
        "Rezende, F. et al. (2021). MAPI: an R package for mapping averaged pairwise "
        "information. *Methods in Ecology and Evolution*, 12, 1305–1310."
    )

# ---------------------------------------------------------------------------
# Assemble document
# ---------------------------------------------------------------------------

lines = [f"# Methods — {project}", ""]
lines.append(
    "> **Note:** Replace `[N_SNPS]` with the SNP count from the pipeline's VCF "
    "statistics output. Replace `[VERSION]` with the installed software version. "
    "All other parameters are filled in from `config/config.yaml`."
)
lines.append("")

for title, body in sections:
    lines.append(f"## {title}")
    lines.append("")
    # Hard-wrap at 90 chars for readability
    for para in body.split("\n\n"):
        lines.append(textwrap.fill(para, width=90))
        lines.append("")

lines.append("## References")
lines.append("")
for ref in refs.values():
    lines.append(f"- {ref}")
lines.append("")

with open(snakemake.output.methods, "w") as fh:
    fh.write("\n".join(lines))
