# Rule for conStruct spatial population structure analysis
# conStruct performs spatial population structure analysis incorporating geography

rule install_construct:
    """
    Install CRAN/r-universe conStruct packages once per workflow run.
    conStruct is not available in the conda channels used here, so the conda
    environment provides compiled dependencies and this rule adds conStruct.
    """
    output:
        touch(".snakemake/construct_installed")
    conda:
        "../envs/construct.yaml"
    shell:
        """
        Rscript --vanilla -e 'lib <- .libPaths()[1]; unlink(list.files(lib, pattern="^00LOCK-", full.names=TRUE), recursive=TRUE, force=TRUE); repos <- c("https://gbradburd.r-universe.dev", "https://cloud.r-project.org"); needed <- c("caroline", "conStruct"); missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing) > 0) install.packages(missing, repos=repos); missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing) > 0) stop("Failed to install required R packages: ", paste(missing, collapse=", "))'
        """


rule construct_analysis:
    """
    Run conStruct analysis for spatial population structure.
    This rule performs the computationally intensive conStruct MCMC analysis
    for different K values (number of layers).
    """
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards),
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_construct.output
    output:
        results_rds = "results/{project}/construct/{project}.construct.K{k}.results.rds",
        layer_proportions = "results/{project}/construct/{project}.construct.K{k}.layer_proportions.txt",
        log_file = "results/{project}/construct/{project}.construct.K{k}.log.txt"
    log:
        "logs/{project}/construct.K{k}.log"
    benchmark:
        "benchmarks/{project}/construct.K{k}.txt"
    params:
        k = lambda wildcards: int(wildcards.k),
        n_chains = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("n_chains", 2),
        n_iterations = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("n_iterations", 50000),
        save_files = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("save.files", True)
    conda:
        "../envs/construct.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["construct"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["construct"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["construct"]["runtime"]
    script:
        "../scripts/construct_analysis.R"


rule construct_choose_k:
    """
    Summarise conStruct MAP log-posterior and layer contributions across K.
    """
    input:
        results_rds = lambda wildcards: expand(
            "results/{project}/construct/{project}.construct.K{k}.results.rds",
            project=wildcards.project,
            k=config["projects"][wildcards.project]["parameters"]["k_values"]
        ),
        install = rules.install_construct.output
    output:
        results_rds = "results/{project}/construct/{project}.construct.chooseK.results.rds",
        choose_k_results = "results/{project}/construct/{project}.construct.chooseK_results.txt",
        lpd_summary = "results/{project}/construct/{project}.construct.lpd_summary.txt",
        layer_contribution_summary = "results/{project}/construct/{project}.construct.layer_contribution_summary.txt"
    log:
        "logs/{project}/construct_choose_k.log"
    benchmark:
        "benchmarks/{project}/construct_choose_k.txt"
    conda:
        "../envs/construct.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["construct"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["construct"]["runtime"]
    script:
        "../scripts/construct_choose_k.R"


rule plot_construct_choose_k:
    """
    Plot conStruct model-comparison scores across K.
    """
    input:
        lpd_summary = rules.construct_choose_k.output.lpd_summary,
        layer_contribution_summary = rules.construct_choose_k.output.layer_contribution_summary
    output:
        lpd_pdf = "results/{project}/construct/plots/{project}.construct.lpd_plot.pdf",
        lpd_rds = "results/{project}/construct/plots/{project}.construct.lpd_plot.rds",
        layer_pdf = "results/{project}/construct/plots/{project}.construct.layer_contribution_plot.pdf",
        layer_rds = "results/{project}/construct/plots/{project}.construct.layer_contribution_plot.rds"
    log:
        "logs/{project}/plot_construct_choose_k.log"
    params:
        width = lambda wildcards: _choose_k_plot_param(wildcards, "width", 10),
        height = lambda wildcards: _choose_k_plot_param(wildcards, "height", 5),
        dpi = lambda wildcards: _choose_k_plot_param(wildcards, "dpi", 300)
    conda:
        "../envs/mapmixture.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_construct_choose_k.R"


rule construct_qmatrix:
    """
    Convert layer proportions to the headerless Q matrix used by mapmixture.
    """
    input:
        layer_proportions = rules.construct_analysis.output.layer_proportions
    output:
        qmatrix = "results/{project}/construct/{project}.construct.K{k}.Qmatrix.txt"
    log:
        "logs/{project}/construct_qmatrix.K{k}.log"
    benchmark:
        "benchmarks/{project}/construct_qmatrix.K{k}.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/construct_qmatrix.R"


rule mapmixture_construct:
    """
    Plot conStruct layer proportions on the shared mapmixture basemap.
    """
    input:
        unpack(_construct_map_inputs),
    output:
        plot = "results/{project}/construct/plots/{project}.construct.K{k}.map.pdf",
        plot_rds = "results/{project}/construct/plots/{project}.construct.K{k}.map.rds"
    wildcard_constraints:
        k = PLOT_K_WILDCARD_CONSTRAINT
    params:
        lambda wildcards: _mapmixture_map_rule_params(wildcards, "construct", "construct")
    log:
        "logs/{project}/mapmixture_construct.K{k}.log"
    benchmark:
        "benchmarks/{project}/mapmixture_construct.K{k}.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"


rule barplot_construct:
    """
    Draw conStruct layer-proportion barplots with the shared mapmixture script.
    """
    input:
        qmatrix = rules.construct_qmatrix.output.qmatrix,
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_mapmixture.output,
    output:
        barplot = "results/{project}/construct/plots/{project}.construct.K{k}.barplot.pdf",
        barplot_rds = "results/{project}/construct/plots/{project}.construct.K{k}.barplot.rds"
    wildcard_constraints:
        k = PLOT_K_WILDCARD_CONSTRAINT
    params:
        lambda wildcards: _barplot_rule_params(wildcards, "construct", "construct")
    log:
        "logs/{project}/construct_barplot.K{k}.log"
    benchmark:
        "benchmarks/{project}/construct_barplot.K{k}.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_barplots.R"

