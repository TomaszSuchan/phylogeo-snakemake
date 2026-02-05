# Rule for conStruct spatial population structure analysis
# conStruct performs spatial population structure analysis incorporating geography

rule construct_analysis:
    """
    Run conStruct analysis for spatial population structure.
    This rule performs the computationally intensive conStruct MCMC analysis
    for different K values (number of layers).
    """
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards),
        indpopdata = rules.generate_popdata.output.indpopdata
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
        n_chains = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("n_chains", 1),
        n_iterations = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("n_iterations", 10000),
        make_freqs = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("make.freqs", True),
        geoDist = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("geoDist", None),
        coords = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("coords", None),
        save_files = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("save.files", True)
    conda:
        "../envs/construct.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["construct"].get("threads", 1)
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["construct"].get("mem_mb", 16000),
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["construct"].get("runtime", 1440)
    script:
        "../scripts/construct_analysis.R"


rule construct_plot:
    """
    Generate conStruct visualization plots from analysis results.
    Creates spatial maps showing layer proportions across geography.
    """
    input:
        results_rds = rules.construct_analysis.output.results_rds,
        indpopdata = rules.generate_popdata.output.indpopdata
    output:
        map_plot = "results/{project}/construct/plots/{project}.construct.K{k}.map.pdf",
        map_plot_rds = "results/{project}/construct/plots/{project}.construct.K{k}.map.rds",
        barplot = "results/{project}/construct/plots/{project}.construct.K{k}.barplot.pdf",
        barplot_rds = "results/{project}/construct/plots/{project}.construct.K{k}.barplot.rds"
    log:
        "logs/{project}/construct_plot.K{k}.log"
    benchmark:
        "benchmarks/{project}/construct_plot.K{k}.txt"
    params:
        crs = lambda wildcards: config["projects"][wildcards.project]["parameters"]["construct"].get("crs", 4326),
        boundary = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("map_boundary", None),
        land_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("land_colour", "#d9d9d9"),
        sea_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("sea_colour", "#deebf7")
    conda:
        "../envs/construct.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"].get("mem_mb", 8000),
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"].get("runtime", 10)
    script:
        "../scripts/construct_plot.R"

