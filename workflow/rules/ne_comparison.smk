"""
Cross-method Ne comparison figure.

Overlays the GONE2 (LD-based, recent generations) and Stairway Plot 2 (folded-SFS,
deep time) trajectories on a shared generation axis, one facet per population, so
the narrow window where the two methods are comparable is visible. The CurrentNe2
point estimate is drawn as a contemporary reference when that analysis is enabled.

One comparison tree per grouping column shared by GONE2 and Stairway Plot 2:
results/{project}/ne_comparison/{grouping}/...
"""


rule ne_comparison_plot:
    input:
        unpack(_ne_comparison_inputs),
    output:
        pdf="results/{project}/ne_comparison/{grouping}/plots/{project}.ne_comparison.pdf",
        rds="results/{project}/ne_comparison/{grouping}/plots/{project}.ne_comparison.rds",
    params:
        pops=lambda wildcards: _ne_comparison_pops(wildcards.project, wildcards.grouping),
        width=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"].get("ne_comparison", {}).get("plot", {}).get("width"),
            45.72,
        ),
        height=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"].get("ne_comparison", {}).get("plot", {}).get("height"),
            22.86,
        ),
        group_sort_by=lambda wildcards: _stairwayplot2_group_setting(
            wildcards.project,
            wildcards.grouping,
            "sort_by",
        ),
    log:
        "logs/{project}/ne_comparison_plot.{grouping}.log"
    benchmark:
        "benchmarks/{project}/ne_comparison_plot_{grouping}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_ne_comparison.R"
