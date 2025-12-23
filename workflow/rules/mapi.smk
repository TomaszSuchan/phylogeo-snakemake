rule mapi_analysis:
    """
    Run MAPI analysis to identify spatial patterns of genetic diversity.
    This rule performs the computationally intensive MAPI calculations including
    permutation tests to generate spatial grid results and significance tails.
    """
    input:
        euclidean_dist="results/{project}/gen_dist/{project}.euclidean_distance.tsv",
        indpopdata="results/{project}/indpopdata.txt",
    output:
        mapi_gpkg="results/{project}/mapi/{project}.mapi_results.gpkg",
        upper_tails_gpkg="results/{project}/mapi/{project}.mapi_upper_tails.gpkg",
        lower_tails_gpkg="results/{project}/mapi/{project}.mapi_lower_tails.gpkg",
    conda:
        "../envs/mapi.yaml"
    log:
        "logs/{project}/mapi_analysis.log",
    benchmark:
        "benchmarks/{project}/mapi_analysis.txt",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["mapi"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["mapi"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["mapi"]["runtime"],
    params:
        n_permutations=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("n_permutations", 1000),
        grid_halfwidth=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("grid_halfwidth", None),
        crs_projected=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("crs_projected", 3857),
        crs_geographic=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("crs_geographic", 4326),
        alpha=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("alpha", 0.05),
    script:
        "../scripts/mapi_analysis.R"


rule mapi_plot:
    """
    Generate MAPI visualization plots from analysis results.
    This lightweight rule creates publication-ready plots using modern ggplot2.
    Can be re-run with different parameters without re-computing MAPI analysis.
    """
    input:
        mapi_gpkg="results/{project}/mapi/{project}.mapi_results.gpkg",
        upper_tails_gpkg="results/{project}/mapi/{project}.mapi_upper_tails.gpkg",
        lower_tails_gpkg="results/{project}/mapi/{project}.mapi_lower_tails.gpkg",
    output:
        mapi_plot="results/{project}/mapi/{project}.mapi_euclidean.pdf",
    conda:
        "../envs/mapi.yaml"
    log:
        "logs/{project}/mapi_plot.log",
    benchmark:
        "benchmarks/{project}/mapi_plot.txt",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    params:
        fill_var=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("fill_var", "avg_value"),
    script:
        "../scripts/mapi_plot.R"
