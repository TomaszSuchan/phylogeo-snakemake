"""
Cross-method Ne comparison figure.

Overlays the GONE2 (LD-based, recent generations) and Stairway Plot 2 (folded-SFS,
deep time) trajectories on a shared generation axis, one facet per population, so
the narrow window where the two methods are comparable is visible. The CurrentNe2
point estimate is drawn as a contemporary reference when that analysis is enabled.
Restricted to populations that both methods estimated: the GONE2 and Stairway Plot 2
population sets are configured independently and need not match.
"""


rule ne_comparison_plot:
    input:
        unpack(_ne_comparison_inputs),
    output:
        pdf="results/{project}/ne_comparison/plots/{project}.ne_comparison.pdf",
        rds="results/{project}/ne_comparison/plots/{project}.ne_comparison.rds",
    params:
        pops=lambda wildcards: _ne_comparison_pops(wildcards.project),
        width=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"].get("ne_comparison", {}).get("plot", {}).get("width"),
            45.72,
        ),
        height=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"].get("ne_comparison", {}).get("plot", {}).get("height"),
            22.86,
        ),
        group_sort_by=lambda wildcards: _easysfs_group_setting(
            wildcards.project,
            config["projects"][wildcards.project]["parameters"]["easysfs"].get("population_column", "Site"),
            "sort_by",
        ),
    log:
        "logs/{project}/ne_comparison_plot.log"
    benchmark:
        "benchmarks/{project}/ne_comparison_plot.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_ne_comparison.R"
