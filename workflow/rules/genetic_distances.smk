"""
Rules for calculating genetic distances between individuals.
These distance matrices can be used for various analyses including MAPI,
clustering, and visualization.
"""

rule kosman_distance:
    """
    Kosman–Leonard (2005) pairwise genetic distance from PLINK dosages (0/1/2).
    Implemented in Python with bed-reader to avoid vcfR/genind OOM on large VCFs.
    """
    input:
        bed=rules.vcf_to_plink_biallelic.output.bed,
        bim=rules.vcf_to_plink_biallelic.output.bim,
        fam=rules.vcf_to_plink_biallelic.output.fam,
    output:
        dist="results/{project}/gen_dist/{project}.kosman_distance.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/{project}/kosman_distance.log",
    benchmark:
        "benchmarks/{project}/kosman_distance.txt",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["runtime"],
    script:
        "../scripts/calculate_kosman_distance.py"


rule euclidean_distance:
    """
    Calculate Euclidean genetic distance from PLINK genotypes.
    Distances are computed from per-SNP dosage values (0, 1, 2 copies
    of the alternate allele) using the existing PLINK files generated
    by preprocessing.
    """
    input:
        bed=rules.vcf_to_plink_biallelic.output.bed,
        bim=rules.vcf_to_plink_biallelic.output.bim,
        fam=rules.vcf_to_plink_biallelic.output.fam,
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


rule p_distance:
    """
    True per-site p-distance from the all-sites (variant + invariant) VCF.
    Computed over biallelic + invariant sites as the pairwise mean of
    |g_i - g_j| / 2 (diploid dosages 0/1/2), where invariant sites contribute 0
    to the numerator and 1 to the denominator, normalising the SNP dissimilarity
    into a true per-site distance. Uses the same reconstructed all-sites VCF as
    pixy, so this analysis requires the ipyrad .loci file.
    """
    input:
        vcf=rules.prepare_invariant_vcf_gz.output.vcf,
    output:
        dist="results/{project}/gen_dist/{project}.p_distance.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/{project}/p_distance.log",
    benchmark:
        "benchmarks/{project}/p_distance.txt",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["gen_dist"]["runtime"],
    script:
        "../scripts/calculate_p_distance.py"


rule average_squared_genetic_difference:
    """
    Calculate the average squared genetic difference matrix.
    Uses the bed2diffs_v2-style formulation:
      similarity = G G^T / n_sites
      diffs_ij = similarity_ii + similarity_jj - 2 * similarity_ij
    """
    input:
        bed=rules.vcf_to_plink_biallelic.output.bed,
        bim=rules.vcf_to_plink_biallelic.output.bim,
        fam=rules.vcf_to_plink_biallelic.output.fam,
    output:
        dist="results/{project}/gen_dist/{project}.avgsquared_distance.tsv",
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

