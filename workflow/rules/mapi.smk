rule mapi_analysis:
    """
    Run MAPI analysis to identify spatial patterns of genetic diversity.
    This rule performs the computationally intensive MAPI calculations including
    permutation tests to generate spatial grid results and significance tails.
    """
    input:
        euclidean_dist=rules.euclidean_distance.output.dist,
        indpopdata=rules.generate_popdata.output.indpopdata,
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
        mapi_gpkg=rules.mapi_analysis.output.mapi_gpkg,
        upper_tails_gpkg=rules.mapi_analysis.output.upper_tails_gpkg,
        lower_tails_gpkg=rules.mapi_analysis.output.lower_tails_gpkg,
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        mapi_plot="results/{project}/mapi/plots/{project}.mapi_euclidean.pdf",
        mapi_plot_rds="results/{project}/mapi/plots/{project}.mapi_euclidean.rds",
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
        width=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("width", 10),
        height=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("height", 8),
        dpi=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("dpi", 300),
        boundary=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("map_boundary", "NULL"),
        crs_plot=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("crs", 4326),
        land_colour=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("land_colour", "#d9d9d9"),
        sea_colour=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("sea_colour", "#deebf7"),
        expand=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("expand", False),
        plot_title=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("plot_title", ""),
        axis_title_size=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("axis_title_size", 10),
        axis_text_size=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("axis_text_size", 8),
        basemap_border=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("basemap_border", True),
        basemap_border_col=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("basemap_border_col", "black"),
        basemap_border_lwd=lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("basemap_border_lwd", 0.1),
        point_size=lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("point_size", 1),
        point_color=lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("point_color", "black"),
        point_alpha=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("point_alpha", 0.6),
        tail_linewidth=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("tail_linewidth", 0.4),
        upper_tail_color=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("upper_tail_color", "#B2182B"),
        lower_tail_color=lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapi"].get("lower_tail_color", "#1B7837"),
    script:
        "../scripts/mapi_plot.R"
