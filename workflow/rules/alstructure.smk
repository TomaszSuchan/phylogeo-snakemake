# ALStructure ancestry estimation, K estimation, and plotting.
#
# ALStructure (Cabreros & Storey 2019) fits the same PSD admixture model as
# STRUCTURE/ADMIXTURE but estimates it by latent-subspace estimation + alternating
# least squares — a spectral, likelihood-free procedure with no MCMC/EM local
# optima, so it converges deterministically and is robust to unbalanced sampling.
# It is the same model family as SCOPE but installs as a pure-R package, and it
# estimates K natively via estimate_d(). Non-spatial: barplots are always produced
# and map plots only when popdata coordinates exist (same gating as ADMIXTURE/sNMF).


def _alstructure_params(wildcards):
    return config["projects"][wildcards.project]["parameters"].get("alstructure", {})


rule install_alstructure:
    """
    Install ALStructure once into the Snakemake conda environment from GitHub.
    Pure-R package (StoreyLab/alstructure); only hard dependency is `svd`,
    provided by the conda env so remotes does not recompile it.
    """
    output:
        touch(".snakemake/alstructure_installed")
    conda:
        "../envs/alstructure.yaml"
    threads: config["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = config["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = config["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        Rscript --vanilla -e 'lib <- .libPaths()[1]; unlink(list.files(lib, pattern="^00LOCK-", full.names=TRUE), recursive=TRUE, force=TRUE); if (!requireNamespace("alstructure", quietly=TRUE, lib.loc=lib)) remotes::install_github("StoreyLab/alstructure", upgrade="never", dependencies=TRUE, build_vignettes=FALSE, lib=lib); if (!requireNamespace("alstructure", quietly=TRUE, lib.loc=lib)) stop("Failed to install required R package: alstructure")'
        """


rule alstructure_genotypes:
    """
    Convert the analysis VCF to an ALStructure-compatible genotype dosage matrix.
    """
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards),
        indpopdata = rules.generate_popdata.output.indpopdata
    output:
        geno_rds = "results/{project}/filtered_data/{project}.alstructure.genotypes.rds"
    log:
        "logs/{project}/alstructure_genotypes.log"
    benchmark:
        "benchmarks/{project}/alstructure_genotypes.txt"
    params:
        ploidy = lambda wildcards: _alstructure_params(wildcards).get("ploidy", 2)
    conda:
        "../envs/alstructure.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["alstructure"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["alstructure"]["runtime"]
    script:
        "../scripts/alstructure_genotypes.R"


rule alstructure_estimate_d:
    """
    Estimate the number of ancestral populations (latent dimension d) and the
    latent-subspace eigenvalue scree across the configured K range.
    """
    input:
        geno_rds = rules.alstructure_genotypes.output.geno_rds,
        install = rules.install_alstructure.output
    output:
        choose_k_results = "results/{project}/alstructure/{project}.alstructure.chooseK_results.txt",
        cv_summary = "results/{project}/alstructure/{project}.alstructure.cv_summary.txt"
    log:
        "logs/{project}/alstructure_estimate_d.log"
    benchmark:
        "benchmarks/{project}/alstructure_estimate_d.txt"
    params:
        k_values = lambda wildcards: config["projects"][wildcards.project]["parameters"]["k_values"]
    conda:
        "../envs/alstructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["alstructure"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["alstructure"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["alstructure"]["runtime"]
    script:
        "../scripts/alstructure_estimate_d.R"


rule plot_alstructure_eigenvalues:
    """
    Plot the ALStructure latent-subspace eigenvalue scree across K (Evanno-style
    ggplot aesthetics); the elbow corroborates the estimate_d value.
    """
    input:
        cv_summary = rules.alstructure_estimate_d.output.cv_summary
    output:
        pdf = "results/{project}/alstructure/plots/{project}.alstructure.cv_plot.pdf",
        rds = "results/{project}/alstructure/plots/{project}.alstructure.cv_plot.rds"
    log:
        "logs/{project}/plot_alstructure_eigenvalues.log"
    params:
        ylab = "Latent-subspace eigenvalue",
        width = lambda wildcards: _choose_k_plot_param(wildcards, "width", 10),
        height = lambda wildcards: _choose_k_plot_param(wildcards, "height", 5),
        dpi = lambda wildcards: _choose_k_plot_param(wildcards, "dpi", 300)
    conda:
        "../envs/alstructure.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_choose_k_score.R"


rule alstructure_analysis:
    """
    Run ALStructure ancestry estimation for a given K.
    """
    input:
        geno_rds = rules.alstructure_genotypes.output.geno_rds,
        install = rules.install_alstructure.output
    output:
        qmatrix = "results/{project}/alstructure/{project}.alstructure.K{k}.Qmatrix.txt",
        results_rds = "results/{project}/alstructure/{project}.alstructure.K{k}.results.rds"
    log:
        "logs/{project}/alstructure.K{k}.log"
    benchmark:
        "benchmarks/{project}/alstructure.K{k}.txt"
    params:
        k = lambda wildcards: int(wildcards.k),
        svd_method = lambda wildcards: _alstructure_params(wildcards).get("svd_method", "base"),
        tolerance = lambda wildcards: _alstructure_params(wildcards).get("tolerance", 1e-05),
        max_iters = lambda wildcards: _alstructure_params(wildcards).get("max_iters", 1000),
        order_method = lambda wildcards: _alstructure_params(wildcards).get("order_method", "ave_admixture"),
        seed = lambda wildcards: _alstructure_params(wildcards).get("seed", 42)
    conda:
        "../envs/alstructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["alstructure"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["alstructure"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["alstructure"]["runtime"]
    script:
        "../scripts/alstructure_analysis.R"


rule mapmixture_alstructure:
    """
    Plot ALStructure ancestry coefficients on the shared mapmixture basemap.
    """
    input:
        unpack(_alstructure_map_inputs),
    output:
        plot = "results/{project}/alstructure/plots/{project}.alstructure.K{k}.map.pdf",
        plot_rds = "results/{project}/alstructure/plots/{project}.alstructure.K{k}.map.rds"
    params:
        lambda wildcards: _mapmixture_map_rule_params(wildcards, "alstructure", "alstructure")
    log:
        "logs/{project}/mapmixture_alstructure.K{k}.log"
    benchmark:
        "benchmarks/{project}/mapmixture_alstructure.K{k}.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"


rule barplot_alstructure:
    """
    Draw ALStructure ancestry barplots with the shared mapmixture barplot script.
    """
    input:
        qmatrix = rules.alstructure_analysis.output.qmatrix,
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_mapmixture.output,
    output:
        barplot = "results/{project}/alstructure/plots/{project}.alstructure.K{k}.barplot.pdf",
        barplot_rds = "results/{project}/alstructure/plots/{project}.alstructure.K{k}.barplot.rds"
    params:
        lambda wildcards: _barplot_rule_params(wildcards, "alstructure", "alstructure")
    log:
        "logs/{project}/alstructure_barplot.K{k}.log"
    benchmark:
        "benchmarks/{project}/alstructure_barplot.K{k}.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_barplots.R"
