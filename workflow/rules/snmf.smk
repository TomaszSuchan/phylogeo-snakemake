# sNMF (LEA) ancestry estimation, choose-K, and plotting.
#
# sNMF is a fast, IBD-robust complement to STRUCTURE/ADMIXTURE: it estimates
# ancestry by sparse non-negative matrix factorisation (Frichot et al. 2014) and
# selects K from the cross-entropy of masked genotypes. It is non-spatial, so
# barplots are always produced and map plots only when popdata coordinates exist
# (same gating as ADMIXTURE/fastStructure).


def _snmf_params(wildcards):
    return config["projects"][wildcards.project]["parameters"].get("snmf", {})


rule snmf_genotypes:
    """
    Convert the analysis VCF to an LEA/sNMF .geno genotype file.
    """
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards),
        indpopdata = rules.generate_popdata.output.indpopdata
    output:
        geno = "results/{project}/filtered_data/{project}.snmf.geno",
        samples = "results/{project}/filtered_data/{project}.snmf.samples.txt"
    log:
        "logs/{project}/snmf_genotypes.log"
    benchmark:
        "benchmarks/{project}/snmf_genotypes.txt"
    params:
        ploidy = lambda wildcards: _snmf_params(wildcards).get("ploidy", 2)
    conda:
        "../envs/snmf.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["snmf"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["snmf"]["runtime"]
    script:
        "../scripts/snmf_genotypes.R"


rule snmf_choose_k:
    """
    Run sNMF across all configured K values and summarise cross-entropy.
    """
    input:
        geno = rules.snmf_genotypes.output.geno
    output:
        results_rds = "results/{project}/snmf/{project}.snmf.chooseK.results.rds",
        choose_k_results = "results/{project}/snmf/{project}.snmf.chooseK_results.txt",
        cv_summary = "results/{project}/snmf/{project}.snmf.cv_summary.txt"
    log:
        "logs/{project}/snmf_choose_k.log"
    benchmark:
        "benchmarks/{project}/snmf_choose_k.txt"
    params:
        k_values = lambda wildcards: config["projects"][wildcards.project]["parameters"]["k_values"],
        repetitions = lambda wildcards: _snmf_params(wildcards).get("repetitions", 10),
        alpha = lambda wildcards: _snmf_params(wildcards).get("alpha", 10),
        tolerance = lambda wildcards: _snmf_params(wildcards).get("tolerance", 1e-05),
        iterations = lambda wildcards: _snmf_params(wildcards).get("iterations", 200),
        ploidy = lambda wildcards: _snmf_params(wildcards).get("ploidy", 2),
        percentage = lambda wildcards: _snmf_params(wildcards).get("percentage", 0.05),
        seed = lambda wildcards: _snmf_params(wildcards).get("seed", 42)
    conda:
        "../envs/snmf.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["snmf"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["snmf"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["snmf"]["runtime"]
    script:
        "../scripts/snmf_choose_k.R"


rule plot_snmf_cv:
    """
    Plot sNMF cross-entropy across K (Evanno-style ggplot aesthetics).
    """
    input:
        cv_summary = rules.snmf_choose_k.output.cv_summary
    output:
        pdf = "results/{project}/snmf/plots/{project}.snmf.cv_plot.pdf",
        rds = "results/{project}/snmf/plots/{project}.snmf.cv_plot.rds"
    log:
        "logs/{project}/plot_snmf_cv.log"
    params:
        ylab = "Cross-entropy",
        crossentropy = True,
        width = lambda wildcards: _choose_k_plot_param(wildcards, "width", 10),
        height = lambda wildcards: _choose_k_plot_param(wildcards, "height", 5),
        dpi = lambda wildcards: _choose_k_plot_param(wildcards, "dpi", 300)
    conda:
        "../envs/snmf.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_choose_k_score.R"


rule snmf_analysis:
    """
    Run sNMF ancestry estimation for a given K (best of several repetitions).
    """
    input:
        geno = rules.snmf_genotypes.output.geno,
        samples = rules.snmf_genotypes.output.samples
    output:
        qmatrix = "results/{project}/snmf/{project}.snmf.K{k}.Qmatrix.txt",
        cross_entropy = "results/{project}/snmf/{project}.snmf.K{k}.cross_entropy.txt",
        results_rds = "results/{project}/snmf/{project}.snmf.K{k}.results.rds"
    log:
        "logs/{project}/snmf.K{k}.log"
    benchmark:
        "benchmarks/{project}/snmf.K{k}.txt"
    params:
        k = lambda wildcards: int(wildcards.k),
        repetitions = lambda wildcards: _snmf_params(wildcards).get("repetitions", 10),
        alpha = lambda wildcards: _snmf_params(wildcards).get("alpha", 10),
        tolerance = lambda wildcards: _snmf_params(wildcards).get("tolerance", 1e-05),
        iterations = lambda wildcards: _snmf_params(wildcards).get("iterations", 200),
        ploidy = lambda wildcards: _snmf_params(wildcards).get("ploidy", 2),
        percentage = lambda wildcards: _snmf_params(wildcards).get("percentage", 0.05),
        seed = lambda wildcards: _snmf_params(wildcards).get("seed", 42)
    conda:
        "../envs/snmf.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["snmf"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["snmf"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["snmf"]["runtime"]
    script:
        "../scripts/snmf_analysis.R"


rule mapmixture_snmf:
    """
    Plot sNMF ancestry coefficients on the shared mapmixture basemap.
    """
    input:
        unpack(_snmf_map_inputs),
    output:
        plot = "results/{project}/snmf/plots/{project}.snmf.K{k}.map.pdf",
        plot_rds = "results/{project}/snmf/plots/{project}.snmf.K{k}.map.rds"
    wildcard_constraints:
        k = PLOT_K_WILDCARD_CONSTRAINT
    params:
        lambda wildcards: _mapmixture_map_rule_params(wildcards, "snmf", "snmf")
    log:
        "logs/{project}/mapmixture_snmf.K{k}.log"
    benchmark:
        "benchmarks/{project}/mapmixture_snmf.K{k}.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"


rule barplot_snmf:
    """
    Draw sNMF ancestry barplots with the shared mapmixture barplot script.
    """
    input:
        qmatrix = "results/{project}/snmf/aligned/{project}.snmf.K{k}.Qmatrix.aligned.txt",
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_mapmixture.output,
    output:
        barplot = "results/{project}/snmf/plots/{project}.snmf.K{k}.barplot.pdf",
        barplot_rds = "results/{project}/snmf/plots/{project}.snmf.K{k}.barplot.rds"
    wildcard_constraints:
        k = PLOT_K_WILDCARD_CONSTRAINT
    params:
        lambda wildcards: _barplot_rule_params(wildcards, "snmf", "snmf")
    log:
        "logs/{project}/snmf_barplot.K{k}.log"
    benchmark:
        "benchmarks/{project}/snmf_barplot.K{k}.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_barplots.R"
