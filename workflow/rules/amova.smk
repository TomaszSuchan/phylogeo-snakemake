"""
Rules for AMOVA (Analysis of Molecular Variance) analysis.
AMOVA partitions genetic variance into hierarchical components based on
configurable stratification levels (e.g., Region/Site).
Missing data is handled by mean imputation.
"""

rule amova:
    """
    Perform AMOVA analysis with mean imputation of missing data.
    Stratification levels are specified in config (e.g., ["Region", "Site"]).
    """
    input:
        vcf=lambda wildcards: get_filtered_vcf_output(wildcards),
        popdata="results/{project}/indpopdata.txt"
    output:
        amova="results/{project}/amova/{project}.amova_results.txt",
        plot="results/{project}/amova/plots/{project}.amova_variance_components.pdf",
        plot_rds="results/{project}/amova/plots/{project}.amova_variance_components.rds"
    params:
        strata=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("amova", {}).get("strata", []),
        nperm=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("amova", {}).get("nperm", 999)
    conda:
        "../envs/r-amova.yaml"
    log:
        "logs/{project}/amova.log"
    benchmark:
        "benchmarks/{project}/amova.txt"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("amova", {}).get("mem_mb", 16000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("amova", {}).get("runtime", 120)
    script:
        "../scripts/perform_amova.R"

