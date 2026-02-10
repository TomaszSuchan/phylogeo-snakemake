"""
Rules for calculating and plotting missingness statistics (individual and locus level).
Missingness is calculated for both:
1. Filtered VCF (before thinning) - {project}.filtered.imiss/lmiss
2. Thinned VCF (after thinning) - {project}.biallelic_snps_thinned.imiss/lmiss

Plots: Histograms for imiss and lmiss (filtered and thinned datasets).

Also calculates VCF statistics: chromosomes, RAD fragments, and variants.
"""

# Missingness from filtered VCF (before thinning)
rule calculate_missing_indv_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        imiss = "results/{project}/stats_vcf/filtered/{project}.filtered.imiss"
    log:
        "logs/{project}/calculate_missing_indv_filtered.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/stats_vcf/filtered/{wildcards.project}.filtered"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --missing-indv \
                 --out {params.out_prefix} &> {log}
        """

rule calculate_missing_loci_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        lmiss = "results/{project}/stats_vcf/filtered/{project}.filtered.lmiss"
    log:
        "logs/{project}/calculate_missing_loci_filtered.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/stats_vcf/filtered/{wildcards.project}.filtered"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --missing-site \
                 --out {params.out_prefix} &> {log}
        """

# Missingness from thinned VCF (after thinning)
rule calculate_missing_indv_thinned:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        imiss = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.imiss"
    log:
        "logs/{project}/calculate_missing_indv_thinned.log"
    params:
        out_prefix = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --missing-indv \
                 --out {params.out_prefix} &> {log}
        """

rule calculate_missing_loci_thinned:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        lmiss = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.lmiss"
    log:
        "logs/{project}/calculate_missing_loci_thinned.log"
    params:
        out_prefix = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --missing-site \
                 --out {params.out_prefix} &> {log}
        """

# Plot imiss histogram - filtered dataset
rule plot_imiss_histogram_filtered:
    input:
        imiss = rules.calculate_missing_indv_filtered.output.imiss
    output:
        pdf = "results/{project}/stats_vcf/filtered/plots/{project}.filtered.imiss_histogram.pdf",
        rds = "results/{project}/stats_vcf/filtered/plots/{project}.filtered.imiss_histogram.rds",
        summary = "results/{project}/stats_vcf/filtered/{project}.filtered.imiss_summary.txt"
    log:
        "logs/{project}/plot_imiss_histogram_filtered.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

# Plot imiss histogram - thinned dataset
rule plot_imiss_histogram_thinned:
    input:
        imiss = rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf = "results/{project}/stats_vcf/thinned/plots/{project}.biallelic_snps.imiss_histogram.pdf",
        rds = "results/{project}/stats_vcf/thinned/plots/{project}.biallelic_snps.imiss_histogram.rds",
        summary = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.imiss_summary.txt"
    log:
        "logs/{project}/plot_imiss_histogram_thinned.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

# Plot lmiss histogram - filtered dataset
rule plot_lmiss_histogram_filtered:
    input:
        lmiss = rules.calculate_missing_loci_filtered.output.lmiss
    output:
        pdf = "results/{project}/stats_vcf/filtered/plots/{project}.filtered.lmiss_histogram.pdf",
        rds = "results/{project}/stats_vcf/filtered/plots/{project}.filtered.lmiss_histogram.rds",
        summary = "results/{project}/stats_vcf/filtered/{project}.filtered.lmiss_summary.txt"
    log:
        "logs/{project}/plot_lmiss_histogram_filtered.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

# Plot lmiss histogram - thinned dataset
rule plot_lmiss_histogram_thinned:
    input:
        lmiss = rules.calculate_missing_loci_thinned.output.lmiss
    output:
        pdf = "results/{project}/stats_vcf/thinned/plots/{project}.biallelic_snps.lmiss_histogram.pdf",
        rds = "results/{project}/stats_vcf/thinned/plots/{project}.biallelic_snps.lmiss_histogram.rds",
        summary = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.lmiss_summary.txt"
    log:
        "logs/{project}/plot_lmiss_histogram_thinned.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

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

