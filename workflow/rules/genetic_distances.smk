"""
Rules for calculating genetic distances between individuals.
These distance matrices can be used for various analyses including MAPI,
clustering, and visualization.
"""

def _gen_dist_resource(wildcards, key):
    resources_cfg = config["projects"][wildcards.project]["parameters"]["resources"]
    gen_dist_cfg = resources_cfg.get("gen_dist", {})
    default_cfg = resources_cfg.get("default", {})
    return gen_dist_cfg.get(key, default_cfg.get(key))


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
    threads: lambda wildcards: _gen_dist_resource(wildcards, "threads")
    resources:
        mem_mb=lambda wildcards: _gen_dist_resource(wildcards, "mem_mb"),
        runtime=lambda wildcards: _gen_dist_resource(wildcards, "runtime"),
    script:
        "../scripts/calculate_euclidean_distance.py"


