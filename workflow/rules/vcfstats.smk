# VCF statistics from original VCF (after user sample subset, before relatedness filtering)
rule calculate_stats_vcf_original:
    input:
        vcf = rules.subset_vcf.output.vcf
    output:
        stats = "results/{project}/stats_vcf/original/{project}.subset.vcf_stats.txt"
    log:
        "logs/{project}/calculate_stats_vcf_original.log"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/calculate_vcf_stats.py"

# VCF statistics from filtered VCF (before thinning)
rule calculate_stats_vcf_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        stats = "results/{project}/stats_vcf/filtered/{project}.filtered.vcf_stats.txt"
    log:
        "logs/{project}/calculate_stats_vcf_filtered.log"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/calculate_vcf_stats.py"

# VCF statistics from thinned VCF (after thinning)
rule calculate_stats_vcf_thinned:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        stats = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.vcf_stats.txt"
    log:
        "logs/{project}/calculate_stats_vcf_thinned.log"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/calculate_vcf_stats.py"

# VCF statistics from original VCF (after user sample subset, before relatedness filtering)
rule chromosome_stats_vcf_original:
    input:
        vcf = rules.subset_vcf.output.vcf
    output:
        stats = "results/{project}/stats_vcf/original/{project}.subset.chromosome_stats.txt"
    log:
        "logs/{project}/chromosome_stats_vcf_original.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        zgrep -v "^#" {input.vcf} \
        | cut -f1 \
        | sort \
        | uniq -c \
        | awk '{{print $2, $1}}' > {output.stats} &> {log}
        """

# VCF statistics from filtered VCF (before thinning)
rule chromosome_stats_vcf_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        stats = "results/{project}/stats_vcf/filtered/{project}.filtered.chromosome_stats.txt"
    log:
        "logs/{project}/chromosome_stats_vcf_filtered.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        zgrep -v "^#" {input.vcf} \
        | cut -f1 \
        | sort \
        | uniq -c \
        | awk '{{print $2, $1}}' > {output.stats} &> {log}
        """

# VCF statistics from thinned VCF (after thinning)
rule chromosome_stats_vcf_thinned:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        stats = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.chromosome_stats.txt"
    log:
        "logs/{project}/chromosome_stats_vcf_thinned.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        zgrep -v "^#" {input.vcf} \
        | cut -f1 \
        | sort \
        | uniq -c \
        | awk '{{print $2, $1}}' > {output.stats} &> {log}
        """

