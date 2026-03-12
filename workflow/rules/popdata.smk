rule generate_popdata:
    input:
        vcf=rules.subset_vcf_after_relatedness.output.vcf
    output:
        indpopdata="results/{project}/indpopdata.txt"
    params:
        popdata=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popdata", ""),
        popmap=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popmap", ""),
        separator=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popseparator", "-")
    log:
        "logs/{project}/generate_popdata.log"
    benchmark:
        "benchmarks/{project}/generate_popdata.txt"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/generate_popdata.py"

rule generate_popdata_all:
    input:
        vcf=rules.subset_vcf.output.vcf
    output:
        indpopdata="results/{project}/indpopdata_all.txt"
    params:
        popdata=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popdata", ""),
        popmap=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popmap", ""),
        separator=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popseparator", "-")
    log:
        "logs/{project}/generate_popdata_all.log"
    benchmark:
        "benchmarks/{project}/generate_popdata_all.txt"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/generate_popdata.py"

rule create_population_summary:
    input:
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        summary="results/{project}/stats_samples/{project}.population_summary.txt"
    log:
        "logs/{project}/create_population_summary.log"
    benchmark:
        "benchmarks/{project}/create_population_summary.txt"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/create_population_summary.py"



# Note: calculate_missing_indv has been moved to missingness.smk
# This rule is kept for backward compatibility with calculate_missing_indv_miss
# but the main rule is now in missingness.smk

rule calculate_missing_indv_miss:
    input:
        vcf = rules.filter_missing_vcf.output.vcf
    output:
        imiss = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps_thinned_miss{miss}.imiss"
    log:
        "logs/{project}/calculate_missing_indv_miss_{miss}.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/stats_vcf/thinned/{wildcards.project}.biallelic_snps_thinned_miss{wildcards.miss}"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        mkdir -p results/{wildcards.project}/stats_vcf/thinned
        vcftools --gzvcf {input.vcf} \
                 --missing-indv \
                 --out {params.out_prefix} &> {log}
        """