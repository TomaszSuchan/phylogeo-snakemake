"""
Generate a manuscript-ready methods.md for one project.

Sources used (all filled in automatically):
  - snakemake.params.analyses   – which tools are enabled
  - snakemake.params.parameters – every config parameter
  - data/mainparams             – STRUCTURE BURNIN / NUMREPS
  - vcf_stats.txt               – variant count, RAD-loci count, sample count
  - depth_summary.txt           – sequencing depth summaries
  - workflow/envs/*.yaml        – software versions

The text is intentionally explicit and explains what each method does, so
that it can be pasted into a manuscript Methods section and trimmed down as
needed. Best-K values are NOT written to the methods (they are results).
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


def parse_depth_summary(path):
    """Return selected values from depth_summary.txt."""
    result = {}
    section = None
    key_map = {
        "Number of individuals": "individuals",
        "Number of sites": "sites",
        "Called genotypes with depth": "called_genotypes",
        "Overall mean depth across called genotypes": "overall_mean_depth",
    }
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line == "Per-individual mean depth":
                    section = "individual"
                    continue
                if line == "Per-site mean depth":
                    section = "site"
                    continue
                if ":" not in line:
                    continue
                label, value = [part.strip() for part in line.split(":", 1)]
                if label in key_map:
                    result[key_map[label]] = value
                elif section and label in {"Mean", "Median", "SD", "Minimum", "Maximum"}:
                    result[f"{section}_{label.lower()}_depth"] = value
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


def k_str(k_list):
    """Return a comma-joined string of sorted K values, e.g. '1, 2, 3, 4, 5'."""
    return ", ".join(str(k) for k in sorted(int(k) for k in k_list))


# ── gather inputs ─────────────────────────────────────────────────────────────

analyses = snakemake.params.analyses        # {tool: bool}
p        = snakemake.params.parameters      # project parameters dict
project  = snakemake.params.project
envs_dir = snakemake.params.envs_dir

mainparams = parse_mainparams(snakemake.input.mainparams)
vcf_stats  = parse_vcf_stats(snakemake.input.vcf_stats)
depth_stats = parse_depth_summary(snakemake.input.depth_summary)
versions   = parse_versions(envs_dir)

n_snps     = vcf_stats.get("variants",      "[N_SNPS]")
n_loci     = vcf_stats.get("rad_fragments", "[N_RAD_LOCI]")
n_samples  = vcf_stats.get("samples",       "[N_SAMPLES]")
overall_depth = depth_stats.get("overall_mean_depth", "[MEAN_DEPTH]")
ind_median_depth = depth_stats.get("individual_median_depth", "[IND_MEDIAN_DEPTH]")
site_median_depth = depth_stats.get("site_median_depth", "[SITE_MEDIAN_DEPTH]")

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
        f"Many model-based clustering and ancestry methods assume that sampled "
        f"individuals are unrelated, because the presence of close relatives can "
        f"create spurious genetic clusters and bias allele-frequency estimates. "
        f"To identify and remove such individuals, pairwise kinship was therefore "
        f"estimated prior to SNP filtering using the KING-robust estimator "
        f"(Manichaikul et al. 2010) as implemented in PLINK 2 version "
        f"{v(versions,'plink2')} (Chang et al. 2015), computed from biallelic "
        f"variants with a minor allele count (MAC) of at least {mac_thr}. "
        f"The KING-robust coefficient estimates kinship directly from genotype "
        f"counts and is robust to population structure. Starting from the most "
        f"connected individuals, samples were greedily excluded until no remaining "
        f"pair had a kinship coefficient exceeding {king_thr}, the conventional "
        f"threshold corresponding to approximately second-degree relatives "
        f"(e.g. half-siblings, avuncular or grandparent-grandchild pairs)."
    )

miss_clause = (f"sites with more than {100*f_missing:.0f}% missing genotype calls "
               f"across individuals were excluded to limit the influence of "
               f"poorly sequenced loci, and " if f_missing < 1.0 else "")
maf_clause  = (f", and a minimum minor allele frequency (MAF) of {maf} was imposed "
               f"to further exclude very rare variants"
               if maf > 0 else "")
filt_parts.append(
    f"SNP genotypes called by ipyrad were then filtered with bcftools version "
    f"{v(versions,'bcftools')} (Danecek et al. 2021). Specifically, {miss_clause}"
    f"variant sites were restricted to biallelic single-nucleotide polymorphisms "
    f"(i.e. sites with exactly two alleles, as required by most downstream "
    f"population-genetic software) with a minor allele count greater than one "
    f"(MAC > 1){maf_clause}. Removing singletons (MAC = 1) discards variants seen "
    f"in only a single allele copy, which carry little information for "
    f"population-level inference and are disproportionately likely to represent "
    f"sequencing or genotyping error."
)

if thinning_strat == "thinning":
    thin_desc = {
        "max_coverage": "the SNP genotyped in the greatest number of individuals "
                        "(maximum sample coverage), with ties broken at random",
        "random":       "a single SNP at random",
        "weighted":     "a single SNP at random, with each SNP weighted by the "
                        "number of individuals in which it was genotyped",
    }.get(thin_meth, thin_meth)
    filt_parts.append(
        f"SNPs located within the same RAD locus are separated by at most the "
        f"length of a sequencing read and share the same underlying genealogy, so "
        f"they are strongly physically linked and are not statistically "
        f"independent. Because methods such as STRUCTURE, ADMIXTURE and PCA assume "
        f"that markers are effectively unlinked, the data were thinned to a single "
        f"SNP per RAD locus by retaining, within each locus, {thin_desc}. "
        f"This produced the **unlinked biallelic SNP dataset** "
        f"({n_snps} SNPs across {n_loci} RAD loci; {n_samples} individuals), "
        f"which was used for all analyses that assume linkage equilibrium."
    )
elif thinning_strat == "ld_pruning":
    filt_parts.append(
        f"Because linked SNPs violate the assumption of marker independence made "
        f"by clustering and ordination methods, the data were pruned for linkage "
        f"disequilibrium (LD) using PLINK version {v(versions,'plink')} "
        f"(Chang et al. 2015). LD pruning was performed in a sliding window of "
        f"{ld_win} SNPs advanced {ld_step} SNPs at a time, recursively removing one "
        f"SNP from each pair whose squared correlation (r2) exceeded {ld_r2} until "
        f"no such pair remained within any window. This produced the "
        f"**LD-pruned biallelic SNP dataset** ({n_snps} SNPs, {n_samples} "
        f"individuals), used for analyses that assume approximate linkage "
        f"equilibrium."
    )
else:
    filt_parts.append(
        f"No thinning or LD pruning was applied; all {n_snps} biallelic SNPs were "
        f"retained as the **all-biallelic SNP dataset** ({n_samples} individuals)."
    )

if analyses.get("pixy", False):
    filt_parts.append(
        f"Estimators of nucleotide diversity and absolute divergence require "
        f"information on invariant (monomorphic) as well as variable sites, because "
        f"excluding invariant sites inflates per-site diversity estimates. A "
        f"separate **all-sites dataset**, containing both variable and invariant "
        f"positions, was therefore reconstructed directly from the ipyrad .loci "
        f"output using a custom Python script. This dataset was deliberately not "
        f"subjected to MAF or MAC filtering, as such filters would remove the "
        f"invariant and rare sites that pixy (Korunes & Samuk 2021) requires for "
        f"unbiased estimation."
    )

filt_parts.append(
    f"Sequencing depth was summarised from the final analysis VCF with vcftools "
    f"version {v(versions,'vcftools')} (Danecek et al. 2011), using --depth to "
    f"estimate mean read depth per individual and --site-mean-depth to estimate "
    f"mean read depth per SNP. In the final {dataset_label}, the mean depth across "
    f"called genotypes was {overall_depth}, with a median individual mean depth of "
    f"{ind_median_depth} and a median site mean depth of {site_median_depth}. "
)

sections.append(("Data filtering and dataset construction",
                 "\n\n".join(filt_parts)))


# 2. Population structure ──────────────────────────────────────────────────────

k_vals  = p.get("k_values", list(range(1, 10)))
k_tested = k_str(k_vals)   # e.g. "1, 2, 3, 4, 5, 6, 7, 8, 9"
struct_parts = []

if analyses.get("structure", False):
    burnin  = mainparams.get("BURNIN",  "100000")
    numreps = mainparams.get("NUMREPS", "200000")
    reps    = p.get("structure", {}).get("replicates", 10)
    struct_parts.append(
        f"Bayesian model-based clustering was performed with STRUCTURE version "
        f"{v(versions,'structure')} (Pritchard et al. 2000) on the {dataset_label} "
        f"({n_snps} SNPs, {n_samples} individuals). STRUCTURE assigns individuals "
        f"probabilistically to a predefined number of genetic clusters (K) by "
        f"seeking cluster allele frequencies and individual membership proportions "
        f"that jointly bring each cluster as close as possible to Hardy-Weinberg "
        f"and linkage equilibrium. Analyses used the admixture model, which allows "
        f"each individual to derive ancestry from more than one cluster "
        f"(NOADMIX = 0), with correlated allele frequencies among clusters "
        f"(FREQSCORR = 1) and the admixture parameter alpha inferred from the data "
        f"(INFERALPHA = 1). For each value of K = {k_tested}, "
        f"{reps} independent replicate runs were performed, each comprising a "
        f"burn-in of {int(burnin):,} Markov chain Monte Carlo (MCMC) iterations "
        f"discarded to allow convergence, followed by {int(numreps):,} iterations "
        f"from which parameters were estimated. Replicate runs at each K were "
        f"aligned to correct for label switching using the FullSearch algorithm in "
        f"the pophelper R package version {v(versions,'r-pophelper')} (Francis 2017). "
        f"The number of clusters best supported by the data was assessed with the "
        f"delta-K method of Evanno et al. (2005), which identifies the K at which "
        f"the second-order rate of change in the log-likelihood is greatest."
    )

if analyses.get("faststructure", False):
    tol   = p.get("faststructure", {}).get("tol",   "10e-6")
    prior = p.get("faststructure", {}).get("prior", "simple")
    struct_parts.append(
        f"Ancestry proportions were estimated with fastStructure version "
        f"{v(versions,'faststructure')} (Raj et al. 2014) on the {dataset_label} "
        f"({n_snps} SNPs, {n_samples} individuals). fastStructure fits the same "
        f"admixture model as STRUCTURE but replaces MCMC sampling with a much "
        f"faster variational Bayesian approximation, making it practical for large "
        f"SNP datasets. The '{prior}' prior on allele frequencies was used, with "
        f"runs iterated to a convergence tolerance of {tol}. The model was fit for "
        f"each value of K = {k_tested}, and the range of K supporting the "
        f"data was determined with the chooseK.py utility distributed with "
        f"fastStructure, which reports the K maximising the marginal likelihood of "
        f"the data and the smallest K that captures the structure in the sample."
    )

if analyses.get("admixture", False):
    struct_parts.append(
        f"Ancestry proportions were estimated by maximum likelihood with ADMIXTURE "
        f"version {v(versions,'admixture')} "
        f"(Alexander et al. 2009) on the {dataset_label} ({n_snps} SNPs, "
        f"{n_samples} individuals). ADMIXTURE uses the same admixture likelihood as "
        f"STRUCTURE but maximises it numerically rather than sampling from the "
        f"posterior, yielding point estimates of ancestry very rapidly. The model "
        f"was fit for each value of K = {k_tested}. For every K, ADMIXTURE's "
        f"built-in five-fold cross-validation procedure was used: the genotype "
        f"matrix is divided into five parts, each in turn masked and predicted from "
        f"the remainder, and the K minimising the resulting cross-validation error "
        f"is taken as the best-supported number of clusters."
    )

if analyses.get("tess3", False):
    tp = p.get("tess3", {})
    reps = tp.get("replicates", 1)
    max_iter = tp.get("max_iteration", 200)
    tol = tp.get("tolerance", 1e-05)
    ploidy = tp.get("ploidy", 2)
    method = tp.get("method", "projected.ls")
    mask = tp.get("mask", 0)
    crossvalid = tp.get("crossvalid", False)
    crossentropy = tp.get("crossentropy", False)
    score_desc = "training RMSE across replicate runs"
    if crossvalid and mask and float(mask) > 0:
        score_desc = (
            "masked cross-validation RMSE"
            if not crossentropy
            else "masked cross-validation cross-entropy"
        )
    elif crossentropy:
        score_desc = "cross-entropy across replicate runs"
    struct_parts.append(
        f"Spatial ancestry coefficients were estimated with tess3r "
        f"(Caye et al. 2016) on the {dataset_label} ({n_snps} SNPs, "
        f"{n_samples} individuals), using sample coordinates from the population "
        f"metadata. Genotypes were converted from VCF GT fields to alternate-allele "
        f"dosage matrices before analysis. tess3r fits a geographically "
        f"constrained non-negative matrix factorisation in which ancestry coefficients "
        f"are regularised toward smooth spatial variation. This makes it useful for "
        f"rapidly visualising spatial ancestry gradients and for comparing spatially "
        f"smoothed ancestry solutions with non-spatial clustering results, while "
        f"recognising that the method smooths Q values rather than explicitly "
        f"modelling relatedness decay under isolation by distance. The model was fit "
        f"for each value of K = {k_tested}, with ploidy = {ploidy}, method "
        f"\"{method}\", {reps} replicate start(s), a maximum of "
        f"{int(max_iter):,} iterations, and a convergence tolerance of {tol}. "
        f"The number of ancestral populations was evaluated from {score_desc} "
        f"across K (tess3r choose-K plot). For each K, ancestry was visualised "
        f"with mapmixture pie maps, tess3r interpolated surfaces, and ggplot maps "
        f"combining ggtess3Q with the shared mapmixture basemap; cross-entropy "
        f"summaries and barplots were also inspected."
    )

if analyses.get("snmf", False):
    sp = p.get("snmf", {})
    reps      = sp.get("repetitions", 10)
    alpha     = sp.get("alpha",       10)
    ploidy    = sp.get("ploidy",      2)
    tol       = sp.get("tolerance",   1e-05)
    iters     = sp.get("iterations",  200)
    pct       = sp.get("percentage",  0.05)
    struct_parts.append(
        f"As a fast complement to the model-based clustering above, ancestry "
        f"proportions were also estimated with sparse non-negative matrix "
        f"factorisation as implemented in sNMF (Frichot et al. 2014) in the LEA R "
        f"package (Frichot & François 2015), on the {dataset_label} ({n_snps} SNPs, "
        f"{n_samples} individuals). Genotypes were converted from VCF GT fields to "
        f"alternate-allele dosages (0–{ploidy}, missing genotypes retained and "
        f"handled natively by sNMF). sNMF estimates individual ancestry coefficients "
        f"and ancestral allele frequencies by regularised least squares rather than "
        f"by fitting an explicit Hardy-Weinberg/linkage-equilibrium likelihood, which "
        f"makes it computationally efficient and comparatively robust to departures "
        f"from these assumptions, including the continuous allele-frequency gradients "
        f"produced by isolation by distance that can lead likelihood-based methods to "
        f"over-split clinal variation. The model was fit for each value of "
        f"K = {k_tested} with ploidy = {ploidy}, a regularisation parameter "
        f"alpha = {alpha}, a convergence tolerance of {tol}, and up to "
        f"{int(iters):,} iterations. For each K, {reps} independent runs were "
        f"performed and the run with the lowest cross-entropy was retained. The "
        f"number of ancestral populations was assessed from the cross-entropy "
        f"criterion, computed on a masked fraction ({pct}) of held-out genotypes "
        f"and minimised across K (Alexander & Lange 2011; Frichot et al. 2014). "
        f"Ancestry coefficients were visualised as mapmixture barplots and, where "
        f"sample coordinates were available, as pie maps on the shared basemap."
    )

if analyses.get("alstructure", False):
    ap = p.get("alstructure", {})
    a_ploidy = ap.get("ploidy", 2)
    a_svd    = ap.get("svd_method", "base")
    a_order  = ap.get("order_method", "ave_admixture")
    struct_parts.append(
        f"Ancestry proportions were additionally estimated with ALStructure "
        f"(Cabreros & Storey 2019), a likelihood-free estimator of the same PSD "
        f"admixture model, on the {dataset_label} ({n_snps} SNPs, {n_samples} "
        f"individuals). Genotypes were converted from VCF GT fields to "
        f"alternate-allele dosages (0–{a_ploidy}), with missing genotypes imputed "
        f"by per-locus means. ALStructure estimates the latent admixture subspace "
        f"from the spectrum of a bias-corrected genotype covariance matrix and then "
        f"recovers ancestral allele frequencies and individual ancestry coefficients "
        f"by alternating least squares; being a deterministic spectral procedure "
        f"rather than MCMC or block-relaxation EM, it avoids random restarts and "
        f"local-optima/convergence issues and is comparatively robust to unbalanced "
        f"sampling. The number of ancestral populations was estimated natively from "
        f"the rank of the latent subspace with estimate_d (Leek 2011), and the "
        f"latent-subspace eigenvalue scree across K = {k_tested} was inspected to "
        f"corroborate it. The model was fit for each K = {k_tested} using the "
        f"\"{a_svd}\" SVD method and \"{a_order}\" cluster ordering. Ancestry "
        f"coefficients were visualised as mapmixture barplots and, where sample "
        f"coordinates were available, as pie maps on the shared basemap."
    )

if analyses.get("evaladmix", False):
    eval_methods = []
    if analyses.get("admixture", False):
        eval_methods.append("ADMIXTURE")
    if analyses.get("faststructure", False):
        eval_methods.append("fastStructure")
    if analyses.get("structure", False):
        eval_methods.append("STRUCTURE")
    if eval_methods:
        methods_str = ", ".join(eval_methods)
        struct_parts.append(
            f"The adequacy of the inferred ancestry models was assessed with "
            f"evalAdmix (Gompert et al. 2014) for {methods_str}. evalAdmix "
            f"compares the residual correlation structure of genotypes not "
            f"explained by the fitted Q and P matrices against that expected "
            f"under the model, helping to identify misspecified K values or "
            f"populations with correlated ancestry not captured by the model."
        )

if analyses.get("dapc", False):
    dp    = p.get("dapc", {})
    n_pca = dp.get("n_pca",     50)
    n_da  = dp.get("n_da",      10)
    # K range mirrors k_values filtered to >= 2 (adegenet requires at least 2 clusters)
    dapc_ks   = [int(k) for k in k_vals if int(k) >= 2]
    k_tested_da = k_str(dapc_ks)
    struct_parts.append(
        f"Discriminant Analysis of Principal Components (DAPC; Jombart et al. 2010), "
        f"implemented in the R package adegenet, was performed on the {dataset_label} "
        f"({n_snps} SNPs, {n_samples} individuals). DAPC is a multivariate approach "
        f"that makes no assumptions about Hardy-Weinberg or linkage equilibrium: it "
        f"first transforms the genotypes with a principal component analysis to remove "
        f"correlations among alleles, then applies a discriminant analysis that "
        f"maximises among-group variation while minimising within-group variation. "
        f"The analysis was run separately for each value of K = {k_tested_da}, "
        f"each time retaining {n_pca} principal components and {n_da} discriminant "
        f"functions. The Bayesian Information Criterion (BIC) computed for each K was "
        f"compared across runs to identify the best-supported number of clusters; "
        f"the number of retained PCA components was kept well below the number of "
        f"individuals to avoid over-fitting."
    )

if analyses.get("spca", False):
    sp = p.get("spca", {})
    nfposi = sp.get("nfposi", 2)
    nfnega = sp.get("nfnega", 2)
    struct_parts.append(
        f"Spatial principal component analysis (sPCA; Jombart et al. 2008), "
        f"implemented in the R package adegenet, was performed on the {dataset_label} "
        f"({n_snps} SNPs, {n_samples} individuals) using sample coordinates from "
        f"the population metadata. sPCA constructs Moran's eigenvector maps from "
        f"geographic positions and finds axes that explain genetic variation while "
        f"separating global structure (positive spatial autocorrelation, broad clines) "
        f"from local structure (negative autocorrelation, fine-scale patches). "
        f"The analysis retained {nfposi} global and {nfnega} local axes. "
        f"Global and local spatial structure was assessed with Monte Carlo tests "
        f"(global.rtest and local.rtest; {sp.get('nperm', 999)} permutations)."
    )

if analyses.get("construct", False):
    cp       = p.get("construct", {})
    n_iter   = cp.get("n_iterations", 10000)
    n_chains = cp.get("n_chains",     1)
    struct_parts.append(
        f"To distinguish discrete population structure from structure arising "
        f"simply through isolation by distance, spatially explicit clustering was "
        f"performed with conStruct (Bradburd et al. 2018) on the {dataset_label} "
        f"({n_snps} SNPs, {n_samples} individuals). This analysis was included "
        f"because standard ancestry-clustering methods can partition continuous "
        f"isolation-by-distance gradients into apparently discrete clusters, often "
        f"assigning individuals near opposite ends of a sampled range to different "
        f"ancestry components and central individuals to apparent admixture. "
        f"Unlike methods that only impose spatial smoothness on ancestry coefficients, "
        f"conStruct models the decay of relatedness with geographic distance within "
        f"each layer and asks how much discrete structure remains after this "
        f"continuous spatial covariance has been accounted for. Models with "
        f"K = {k_tested} spatial layers were fit, each using {n_chains} MCMC "
        f"chain(s) of {n_iter:,} iterations. Spatial and non-spatial model fits, "
        f"cross-validation scores, and the contribution of each layer to the "
        f"explained covariance were compared to identify layer counts that improved "
        f"prediction without over-interpreting clinal variation as discrete "
        f"population structure."
    )

if analyses.get("treemix", False):
    tm = p.get("treemix", {})
    pop_col = tm.get("population_column", "Site")
    migration_edges = tm.get("migration_edges", [0])
    if isinstance(migration_edges, int):
        migration_edges = [migration_edges]
    migration_str = ", ".join(str(int(m)) for m in migration_edges)
    optm_cfg = tm.get("optm", {})
    optm_method = optm_cfg.get("method", "Evanno")
    optm_reps = optm_cfg.get("replicates", 1)
    root_value = tm.get("root", None)
    root_clause = (
        f"The starting tree was rooted with '{root_value}' before graph search. "
        if root_value not in [None, "", "null", "NULL"]
        else "No explicit root was configured, so graph searches were run without the recommended outgroup root. "
    )
    mlno_value = tm.get("mlno", True)
    allmigs_value = tm.get("allmigs", False)
    bootstrap_cfg = tm.get("bootstrap", {})
    if isinstance(bootstrap_cfg, bool):
        bootstrap_cfg = {"enabled": bootstrap_cfg, "migration_edges": [], "replicates": 100}
    bootstrap_enabled = bool(bootstrap_cfg.get("enabled", False))
    bootstrap_edges = bootstrap_cfg.get("migration_edges", [])
    if isinstance(bootstrap_edges, int):
        bootstrap_edges = [bootstrap_edges]
    bootstrap_edges_str = ", ".join(str(edge) for edge in bootstrap_edges) if bootstrap_edges else "none"
    bootstrap_reps = int(bootstrap_cfg.get("replicates", 100))

    def _orientagraph_flag_text(flag, value):
        if value in [None, False, "", "null", "NULL"]:
            return f"{flag} was not used"
        if value is True:
            return flag
        return f"{flag} {value}"

    struct_parts.append(
        f"Historical relationships among populations were inferred with "
        f"OrientAGraph version {v(versions,'orientagraph')} (Molloy et al. 2021), "
        f"which augments the TreeMix framework (Pickrell & Pritchard 2012), from "
        f"the {dataset_label} ({n_snps} SNPs, {n_samples} individuals). Individuals "
        f"were grouped into populations using the '{pop_col}' column of "
        f"indpopdata, and reference and alternate alleles were counted within "
        f"populations at each SNP to produce TreeMix-compatible input. "
        f"Graphs were fit with migration-edge counts m = {migration_str}, using "
        f"a block size of {tm.get('k', tm.get('block_size', 1000))} SNPs to account for linkage "
        f"when estimating the covariance matrix. {root_clause}"
        f"Maximum likelihood network orientation (MLNO) was configured as "
        f"`{_orientagraph_flag_text('-mlno', mlno_value)}`, and exhaustive "
        f"migration-edge search was configured as "
        f"`{_orientagraph_flag_text('-allmigs', allmigs_value)}`. For each migration-edge count, "
        f"bootstrap resampling was "
        f"{f'enabled for selected migration-edge counts m = {bootstrap_edges_str} with {bootstrap_reps} independent replicate(s) per m' if bootstrap_enabled else 'not enabled for the default maximum-likelihood graph fit'}. "
        f"The inferred graph and observed-minus-model covariance residuals were "
        f"plotted. Migration-edge support was summarized with OptM version "
        f"{v(versions,'r-optm')} using the {optm_method} method across "
        f"{optm_reps} independent run(s) per migration-edge count. OptM estimates "
        f"the most useful value of m from the second-order rate of change in "
        f"log-likelihood (delta-m), which prioritizes migration edges that improve "
        f"model fit rather than proving the true number of historical migration "
        f"events. Final log likelihoods were also compared across migration-edge "
        f"counts."
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
    SVD_NAMES = {0: "implicitly restarted Arnoldi (IRAM)",
                 1: "single-pass randomised SVD",
                 2: "window-based randomised SVD",
                 3: "full singular value decomposition"}
    svd_desc = SVD_NAMES.get(int(svd_m), f"method {svd_m}")
    pca_parts.append(
        f"Principal component analysis (PCA) was used as a model-free summary of "
        f"genetic variation, projecting individuals onto orthogonal axes that "
        f"successively capture the largest fractions of allele-frequency variance "
        f"and typically separating individuals by ancestry along the leading axes. "
        f"PCA was performed with PCAone version {v(versions,'pcaone')} "
        f"(Zhang & Meisner 2023) on the {dataset_label} ({n_snps} SNPs, "
        f"{n_samples} individuals), extracting the {n_pc} leading principal "
        f"components via {svd_desc}."
    )

if analyses.get("pcaone_emu", False) or analyses.get("emupca", False):
    pp   = p.get("PCAone", {})
    n_pc = pp.get("PCnum", 10)
    pca_parts.append(
        f"Because missing genotypes can distort principal components, PCA was also "
        f"run in EMU mode in PCAone version {v(versions,'pcaone')}, which uses an "
        f"iterative expectation-maximisation algorithm to account for missing data "
        f"while estimating the {n_pc} leading principal components, and the results "
        f"were compared with the standard PCA above."
    )

if analyses.get("pcoa", False):
    pca_parts.append(
        f"Principal coordinates analysis (PCoA), the distance-based analogue of "
        f"PCA, was performed on the pairwise genetic distance matrices described "
        f"below using the cmdscale function in R, embedding individuals in a "
        f"low-dimensional space that preserves their pairwise genetic distances as "
        f"faithfully as possible."
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
    STAT_LABELS = {"pi":  "nucleotide diversity (pi, the average number of "
                          "pairwise differences per site within populations)",
                   "fst": "Fst (the proportion of total genetic variance due to "
                          "differences between populations)",
                   "dxy": "Dxy (absolute nucleotide divergence between populations)"}
    stats_str = ", ".join(STAT_LABELS.get(s, s) for s in stats)
    verb = "was" if len(stats) == 1 else "were"
    div_parts.append(
        f"Within- and between-population genetic diversity was quantified with "
        f"pixy version {v(versions,'pixy')} (Korunes & Samuk 2021), which estimates "
        f"{stats_str}. pixy computes these statistics directly from the all-sites "
        f"dataset, counting both variable and invariant positions so that missing "
        f"data do not bias per-site estimates. Statistics {verb} calculated in "
        f"non-overlapping windows of {win:,} bp, and population-level point "
        f"estimates with 95% confidence intervals were obtained by bootstrapping "
        f"across windows ({nboot:,} replicates)."
    )

if analyses.get("genome_scan", False):
    gs     = p.get("genome_scan", {})
    pop1   = gs.get("pop1",  "")
    pop2   = gs.get("pop2",  "")
    w_fst  = gs.get("window_size_fst",  100)
    w_pi   = gs.get("window_size_pi",  1_000_000)
    w_dxy  = gs.get("window_size_dxy", 1_000_000)
    div_parts.append(
        f"Genome-wide scans of differentiation between {pop1} and {pop2} were "
        f"performed with pixy version {v(versions,'pixy')} on an invariant-site "
        f"VCF reconstructed from the original ipyrad loci file, retaining only "
        f"individuals assigned to the two focal groups. Sliding-window "
        f"estimates of F_ST (window size {w_fst:,} bp), nucleotide diversity "
        f"pi (window size {w_pi:,} bp), and absolute divergence Dxy "
        f"(window size {w_dxy:,} bp) were plotted along the genome to "
        f"localise regions of elevated differentiation or diversity."
    )

if analyses.get("amova", False):
    ap     = p.get("amova", {})
    strata = ap.get("strata", [])
    nperm  = ap.get("nperm",  999)
    strata_str = " and ".join(strata) if strata else "[hierarchical groupings]"
    div_parts.append(
        f"The hierarchical partitioning of genetic variance was tested with an "
        f"analysis of molecular variance (AMOVA; Excoffier et al. 1992), which "
        f"decomposes the total genetic variance into components attributable to "
        f"different levels of a sampling hierarchy. AMOVA was carried out in R with "
        f"the poppr package (Kamvar et al. 2014) on the {dataset_label}, "
        f"partitioning variance among {strata_str}. The statistical significance of "
        f"each variance component was assessed with {nperm} random permutations of "
        f"individuals among groups."
    )

if div_parts:
    sections.append(("Genetic diversity and differentiation",
                     "\n\n".join(div_parts)))


# 5. Genetic distances and networks ───────────────────────────────────────────

dist_parts = []

if analyses.get("gen_dist", False):
    dist_parts.append(
        f"To summarise overall genetic similarity among individuals independently "
        f"of any population model, four complementary pairwise genetic distance "
        f"matrices were computed from the {dataset_label} ({n_snps} SNPs, "
        f"{n_samples} individuals) using custom Python scripts that read PLINK "
        f"binary genotype files via the bed-reader library: "
        f"(i) the Kosman-Leonard distance (Kosman & Leonard 2005), a distance "
        f"designed for diploid codominant markers that scores heterozygotes as "
        f"intermediate between the two homozygotes; "
        f"(ii) the Euclidean distance between per-SNP allele-dosage vectors "
        f"(coded 0, 1, 2 for the number of copies of the alternate allele), with "
        f"missing genotypes replaced by the per-SNP mean; "
        f"(iii) the p-distance, the proportion of allelic differences between two "
        f"individuals averaged over all loci scored in both; and "
        f"(iv) the average squared genotype difference (the bed2diffs formulation "
        f"used by EEMS), which underlies several landscape-genetic methods."
    )

if analyses.get("neighbornet", False):
    nn_color  = p.get("neighbornet", {}).get("color_by", ["Site"])
    color_str = ", ".join(nn_color) if isinstance(nn_color, list) else nn_color
    # If the genetic-distance section above already defined the p-distance, refer
    # back to it; otherwise describe its computation here so this paragraph is
    # self-contained.
    if analyses.get("gen_dist", False):
        pdist_clause = "the p-distance matrix described above"
    else:
        pdist_clause = (
            f"a pairwise p-distance matrix, computed from the {dataset_label} "
            f"({n_snps} SNPs, {n_samples} individuals) with a custom Python script "
            f"reading PLINK binary genotype files via the bed-reader library as the "
            f"proportion of allelic differences between two individuals averaged "
            f"over all loci genotyped in both"
        )
    dist_parts.append(
        f"To visualise relationships among individuals while allowing for "
        f"conflicting or reticulate signal (which a strictly bifurcating tree "
        f"cannot represent), a NeighborNet split network (Bryant & Moulton 2004) "
        f"was constructed with the phangorn package (Schliep 2011) in R and drawn "
        f"with the tanggle package, with tips coloured by {color_str}. The network "
        f"was built from {pdist_clause}."
    )

if dist_parts:
    sections.append(("Genetic distances and networks",
                     "\n\n".join(dist_parts)))


# 6. Phylogenetics ────────────────────────────────────────────────────────────

phy_parts = []

for iq_key, label in [("iqtree",        ""),
                      ("iqtree_trimal", " with trimAl alignment trimming"),
                      ("iqtree_robust", " (robust-phylogeny mode)")]:
    if not analyses.get(iq_key, False):
        continue
    ip       = p.get("iqtree", {})
    model    = ip.get("model",             "MFP")
    boots    = ip.get("bootstraps",        1000)
    outgroup = ip.get("outgroup",          None)
    supp_thr = ip.get("support_threshold", 70)

    model_clause = (
        "the best-fit nucleotide substitution model selected automatically by "
        "ModelFinder (Kalyaanamoorthy et al. 2017) according to the Bayesian "
        "Information Criterion"
        if model == "MFP" else f"the {model} substitution model"
    )
    root_clause = (f"rooted with {outgroup} as the outgroup"
                   if outgroup else "midpoint-rooted")

    trimal_clause = ""
    if iq_key == "iqtree_trimal":
        gt = p.get("trimal", {}).get("gt", 0.95)
        trimal_clause = (
            f" Prior to tree inference, poorly aligned and gap-rich columns were "
            f"removed with trimAl version {v(versions,'trimal')} "
            f"(Capella-Gutierrez et al. 2009) using a gap threshold of {gt}, which "
            f"discards alignment columns present in fewer than {100*gt:.0f}% of "
            f"individuals."
        )

    phy_parts.append(
        f"A maximum-likelihood phylogeny{label} was inferred from the concatenated "
        f"RAD-seq sequence alignment exported by ipyrad using IQ-TREE 2 version "
        f"{v(versions,'iqtree')} (Minh et al. 2020) under {model_clause}.{trimal_clause} "
        f"The tree maximising the likelihood of the alignment under this model was "
        f"retained, and branch support was estimated from {boots:,} ultrafast "
        f"bootstrap replicates (Hoang et al. 2018), a resampling procedure that "
        f"quantifies how consistently each clade is recovered; nodes with less than "
        f"{supp_thr}% support are not shown. The tree was {root_clause} for display."
    )

if analyses.get("fineradstructure", False):
    fr  = p.get("fineradstructure", {})
    cl  = fr.get("cluster", {})
    mc  = cl.get("mcmc_iterations", 1000000)
    bi  = cl.get("burnin",          1000000)
    th  = cl.get("thinning",        1000)
    mct = fr.get("tree", {}).get("mcmc_iterations", 10000)
    frs_ver = v(versions, "fineradstructure")
    # f_missing is the only variant filter that touches the fineRADstructure input
    # (a site-level missingness filter; it removes gappy loci, not rare alleles).
    frs_missing_clause = (
        f"after removing sites with more than {100*f_missing:.0f}% missing "
        f"genotype calls"
        if f_missing < 1.0
        else "without any missingness filtering"
    )
    phy_parts.append(
        f"Fine-scale population structure was additionally inferred from shared "
        f"co-ancestry with fineRADstructure version {frs_ver} (Malinsky et al. "
        f"2018), a method tailored to RAD-seq data that exploits the haplotype "
        f"(phase) information available within each RAD locus and is sensitive to "
        f"recent shared ancestry. Because this co-ancestry signal is carried "
        f"disproportionately by rare alleles, fineRADstructure was run on the "
        f"variant set obtained {frs_missing_clause}, retaining all remaining SNPs "
        f"(including singletons and multiallelic sites) without the minor allele "
        f"count, minor allele frequency, or one-SNP-per-locus thinning filters "
        f"applied to the other analyses, since those filters would remove the "
        f"variants most informative for recent co-ancestry. "
        f"RADpainter first summarised, for every pair of "
        f"individuals, the extent to which they share their most closely related "
        f"haplotypes, producing a pairwise co-ancestry matrix ({n_samples} "
        f"individuals). fineSTRUCTURE (Lawson et al. 2012) then clustered "
        f"individuals from this matrix using MCMC ({mc:,} sampling iterations "
        f"following a {bi:,}-iteration burn-in, sampled every {th:,} iterations), "
        f"and a co-ancestry tree relating the inferred clusters was built from "
        f"{mct:,} tree-building iterations."
    )

if phy_parts:
    sections.append(("Phylogenetic and co-ancestry analysis",
                     "\n\n".join(phy_parts)))


# 7. Relatedness and ROH ──────────────────────────────────────────────────────

rel_parts = []

if analyses.get("relatedness", False):
    rel_parts.append(
        f"Pairwise relatedness among the {n_samples} individuals was estimated "
        f"from the {dataset_label} using four complementary estimators, so that "
        f"conclusions did not depend on the assumptions of any single method: "
        f"(i) the Ajk genomic relatedness statistic of Yang et al. (2010), computed "
        f"with vcftools version {v(versions,'vcftools')} (Danecek et al. 2011) "
        f"using --relatedness; "
        f"(ii) the kinship coefficient of Manichaikul et al. (2010), which is more "
        f"robust to population structure, computed with vcftools --relatedness2; "
        f"(iii) the identity-by-descent proportion (PI_HAT) from PLINK version "
        f"{v(versions,'plink')} --genome; and "
        f"(iv) the KING-robust kinship coefficient from PLINK 2 "
        f"--make-king-table (Manichaikul et al. 2010). Together these statistics "
        f"identify duplicate samples, close relatives, and broader patterns of "
        f"background relatedness."
    )

if analyses.get("roh", False):
    roh_group = p.get("roh", {}).get("group_by", ["Site"])
    group_str = " and ".join(roh_group) if isinstance(roh_group, list) else roh_group
    rel_parts.append(
        f"Runs of homozygosity (ROH) - long stretches of consecutive homozygous "
        f"genotypes that arise when an individual inherits identical haplotypes "
        f"from both parents and therefore signal recent inbreeding - were detected "
        f"with the roh model of bcftools version {v(versions,'bcftools')} "
        f"(Danecek et al. 2021). Because ROH detection relies on the density of "
        f"consecutive genotypes along each locus, ROH were called from the variant "
        f"set after sample subsetting, relatedness filtering, and any per-site "
        f"missingness filtering, but before the biallelic, minor-allele-count, and "
        f"one-SNP-per-locus thinning filters applied to the other analyses, since "
        f"thinning to one SNP per RAD locus would remove the consecutive sites that "
        f"ROH calling requires. The number and total length "
        f"of ROH were summarised per individual and compared across {group_str} "
        f"groups, and the proportion of the genome falling within ROH (F_ROH) was "
        f"used as a genomic estimate of the individual inbreeding coefficient."
    )

if rel_parts:
    sections.append(("Relatedness and inbreeding",
                     "\n\n".join(rel_parts)))


# 8. Spatial genetics ─────────────────────────────────────────────────────────

if analyses.get("mapi", False):
    mp        = p.get("mapi", {})
    n_perm    = mp.get("n_permutations", 1000)
    halfwidth = mp.get("grid_halfwidth")
    alpha     = mp.get("alpha",           0.05)
    crs_proj  = mp.get("crs_projected",   "EPSG:8857")
    mapi_ver  = v(versions, "r-mapi")
    if halfwidth is None:
        grid_desc = (
            "a hexagonal grid with cell half-width auto-estimated from sample "
            "coordinates (MAPI_GridAuto)"
        )
    else:
        grid_desc = f"a hexagonal grid of cell half-width {halfwidth:,} m"
    sections.append((
        "Spatial genetic structure",
        f"Geographic patterns in genetic differentiation were mapped with MAPI "
        f"(Mapping Averaged Pairwise Information; Piry et al. 2016) as implemented "
        f"in the R package mapi version {mapi_ver}. MAPI overlays a hexagonal grid "
        f"on the study area and, for every pair of individuals, attributes their "
        f"pairwise genetic distance to all grid cells intersected by an ellipse "
        f"joining their locations; averaging these contributions per cell yields a "
        f"continuous surface highlighting areas of unusually high differentiation "
        f"(barriers to gene flow) or unusually low differentiation (corridors). "
        f"The Euclidean genetic distance matrix ({n_samples} individuals) was used "
        f"as input, with sample coordinates projected to {crs_proj} and {grid_desc}. "
        f"The significance of the upper and lower tails of the cell-value "
        f"distribution was assessed by {n_perm:,} permutations of the genetic "
        f"distances among samples (significance level alpha = {alpha})."
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

if analyses.get("tess3", False):
    refs["tess3"] = (
        "Caye, K., Deist, T.M., Martins, H., Michel, O. & François, O. (2016). "
        "TESS3: fast inference of spatial population structure and genome scans "
        "for selection. *Molecular Ecology Resources*, 16, 540–548."
    )

if analyses.get("snmf", False):
    refs["snmf"] = (
        "Frichot, E., Mathieu, F., Trouillon, T., Bouchard, G. & François, O. "
        "(2014). Fast and efficient estimation of individual ancestry "
        "coefficients. *Genetics*, 196, 973–983."
    )
    refs["lea"] = (
        "Frichot, E. & François, O. (2015). LEA: an R package for landscape and "
        "ecological association studies. *Methods in Ecology and Evolution*, 6, "
        "925–929."
    )
    refs["crossentropy"] = (
        "Alexander, D.H. & Lange, K. (2011). Enhancements to the ADMIXTURE "
        "algorithm for individual ancestry estimation. *BMC Bioinformatics*, "
        "12, 246."
    )

if analyses.get("alstructure", False):
    refs["alstructure"] = (
        "Cabreros, I. & Storey, J.D. (2019). A likelihood-free estimator of "
        "population structure bridging admixture models and principal components "
        "analysis. *Genetics*, 212, 1009–1029."
    )
    refs["estimate_d"] = (
        "Leek, J.T. (2011). Asymptotic conditional singular value decomposition "
        "for high-dimensional genomic data. *Biometrics*, 67, 344–352."
    )

if analyses.get("evaladmix", False):
    refs["evaladmix"] = (
        "Gompert, Z., Buerkle, C.A. & Parchman, T.L. (2014). Genomic evidence "
        "for speciation and gene flow between Timema cristinae stick insects. "
        "*Molecular Ecology*, 23, 1483–1499."
    )

if analyses.get("dapc", False):
    refs["dapc"] = (
        "Jombart, T., Devillard, S. & Balloux, F. (2010). Discriminant analysis of "
        "principal components: a new method for the analysis of genetically structured "
        "populations. *BMC Genetics*, 11, 94."
    )

if analyses.get("spca", False):
    refs["spca"] = (
        "Jombart, T., Devillard, S., Dufour, A.B. & Pontier, D. (2008). Revealing "
        "cryptic spatial patterns in genetic variability by a new multivariate method. "
        "*BMC Genetics*, 9, 461."
    )

if analyses.get("construct", False):
    refs["construct"] = (
        "Bradburd, G.S., Coop, G.M. & Ralph, P.L. (2018). Inferring continuous and "
        "discrete population genetic structure across space. *Genetics*, 210, 33–52."
    )

if analyses.get("treemix", False):
    refs["orientagraph"] = (
        "Molloy, E.K., Durvasula, A. & Sankararaman, S. (2021). Advancing "
        "admixture graph estimation via maximum likelihood network orientation. "
        "*Bioinformatics*, 37(Supplement_1), i142–i150. "
        "https://doi.org/10.1093/bioinformatics/btab267"
    )
    refs["treemix"] = (
        "Pickrell, J.K. & Pritchard, J.K. (2012). Inference of population splits "
        "and mixtures from genome-wide allele frequency data. *PLOS Genetics*, "
        "8, e1002967."
    )
    refs["optm"] = (
        "Fitak, R.R. (2021). OptM: estimating the optimal number of migration edges "
        "on population trees using Treemix. *Biology Methods and Protocols*, "
        "6(1), bpab017. https://doi.org/10.1093/biomethods/bpab017"
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
        "Capella-Gutiérrez, S., Silla-Martínez, J.M. & Gabaldón, T. (2009). "
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
        "Piry, S. et al. (2016). Mapping Averaged Pairwise Information (MAPI): a new "
        "exploratory tool to uncover spatial structure. *Methods in Ecology and "
        "Evolution*, 7, 1463–1475."
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
