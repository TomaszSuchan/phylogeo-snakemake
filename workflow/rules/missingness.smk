"""
Rules for calculating and plotting missingness statistics (individual and locus level).
Missingness is calculated for both:
1. Filtered VCF (before thinning) - {project}.filtered.imiss/lmiss
2. Thinned VCF (after thinning) - {project}.biallelic_snps_thinned.imiss/lmiss

Plots: Histograms for imiss and lmiss (filtered and thinned datasets).
"""

# Missingness from filtered VCF (before thinning)
rule calculate_missing_indv_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        imiss = "results/{project}/missingness_data/filtered/{project}.filtered.imiss"
    log:
        "logs/{project}/calculate_missing_indv_filtered.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/missingness_data/filtered/{wildcards.project}.filtered"
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
        lmiss = "results/{project}/missingness_data/filtered/{project}.filtered.lmiss"
    log:
        "logs/{project}/calculate_missing_loci_filtered.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/missingness_data/filtered/{wildcards.project}.filtered"
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
        vcf = rules.thin_vcf.output.vcf
    output:
        imiss = "results/{project}/missingness_data/thinned/{project}.biallelic_snps_thinned.imiss"
    log:
        "logs/{project}/calculate_missing_indv_thinned.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/missingness_data/thinned/{wildcards.project}.biallelic_snps_thinned"
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
        vcf = rules.thin_vcf.output.vcf
    output:
        lmiss = "results/{project}/missingness_data/thinned/{project}.biallelic_snps_thinned.lmiss"
    log:
        "logs/{project}/calculate_missing_loci_thinned.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/missingness_data/thinned/{wildcards.project}.biallelic_snps_thinned"
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
        pdf = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_histogram.pdf",
        rds = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_histogram.rds"
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
        pdf = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_histogram.pdf",
        rds = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_histogram.rds"
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
        pdf = "results/{project}/missingness_data/filtered/plots/{project}.filtered.lmiss_histogram.pdf",
        rds = "results/{project}/missingness_data/filtered/plots/{project}.filtered.lmiss_histogram.rds"
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
        pdf = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.lmiss_histogram.pdf",
        rds = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.lmiss_histogram.rds"
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

