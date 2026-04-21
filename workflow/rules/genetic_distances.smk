"""
Rules for calculating genetic distances between individuals.
These distance matrices can be used for various analyses including MAPI,
clustering, and visualization.
"""

rule euclidean_distance:
    """
    Calculate Euclidean genetic distance from PLINK genotypes.
    Distances are computed from per-SNP dosage values (0, 1, 2 copies
    of the alternate allele) using the existing PLINK files generated
    by preprocessing.
    """
    input:
        bed=rules.vcf_to_plink.output.bed,
        bim=rules.vcf_to_plink.output.bim,
        fam=rules.vcf_to_plink.output.fam,
    output:
        dist="results/{project}/gen_dist/{project}.euclidean_distance.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/{project}/euclidean_distance.log",
    benchmark:
        "benchmarks/{project}/euclidean_distance.txt",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["runtime"],
    script:
        "../scripts/calculate_euclidean_distance.py"


rule average_squared_genetic_difference:
    """
    Calculate the average squared genetic difference matrix.
    Uses the bed2diffs_v2-style formulation:
      similarity = G G^T / n_sites
      diffs_ij = similarity_ii + similarity_jj - 2 * similarity_ij
    """
    input:
        bed=rules.vcf_to_plink.output.bed,
        bim=rules.vcf_to_plink.output.bim,
        fam=rules.vcf_to_plink.output.fam,
    output:
        dist="results/{project}/gen_dist/{project}.average_squared_genetic_difference.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/{project}/average_squared_genetic_difference.log",
    benchmark:
        "benchmarks/{project}/average_squared_genetic_difference.txt",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["runtime"],
    script:
        "../scripts/calculate_average_squared_genetic_difference.py"

