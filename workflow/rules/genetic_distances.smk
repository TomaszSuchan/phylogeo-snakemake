"""
Rules for calculating genetic distances between individuals.
These distance matrices can be used for various analyses including MAPI,
clustering, and visualization.
"""

rule euclidean_distance:
    """
    Calculate Euclidean genetic distance from VCF.
    Euclidean distance is calculated from the genotype matrix after
    converting to allele counts (0, 1, 2).
    """
    input:
        vcf=lambda wildcards: get_filtered_vcf_output(wildcards),
    output:
        dist="results/{project}/gen_dist/{project}.euclidean_distance.tsv",
    conda:
        "../envs/r-popgenreport.yaml"
    log:
        "logs/{project}/euclidean_distance.log",
    benchmark:
        "benchmarks/{project}/euclidean_distance.txt",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/calculate_euclidean_distance.R"
