# Quality Control Rules
# This file contains rules for QC analyses of VCF data

# Rule to generate comprehensive VCF assembly statistics
rule qc_assembly:
    input:
        vcf="results/{project}/filtered_data/{project}.raw_sorted.vcf.gz"
    output:
        report="results/{project}/qc/{project}.assembly_stats.txt"
    log:
        "logs/{project}/vcf_assembly_stats.log"
    benchmark:
        "benchmarks/{project}/vcf_assembly_stats.txt"
    params:
        id_pattern=lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("id_pattern", r"loc(\d+)_")
    conda:
        "../envs/vcfpy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/vcf_assembly_stats.py \
            --vcf {input.vcf} \
            --out {output.report} \
            --id-pattern '{params.id_pattern}' &> {log}
        """


# Rule to run whoa (Where's my Heterozygotes at?) QC analysis
# Evaluates genotyping accuracy by examining heterozygote miscall rates
rule qc_whoa:
    input:
        vcf="results/{project}/filtered_data/{project}.raw_sorted.vcf.gz"
    output:
        posteriors="results/{project}/qc/{project}.whoa_posteriors.rds",
        plot="results/{project}/qc/{project}.whoa_plot.pdf",
        report="results/{project}/qc/{project}.whoa_report.txt"
    log:
        "logs/{project}/whoa_qc.log"
    benchmark:
        "benchmarks/{project}/whoa_qc.txt"
    params:
        min_bin=lambda wildcards: config["projects"][wildcards.project]["parameters"]["whoa"].get("min_bin", 1000)
    conda:
        "../envs/r-whoa.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        Rscript workflow/scripts/run_whoa_qc.R \
            --vcf {input.vcf} \
            --posteriors {output.posteriors} \
            --plot {output.plot} \
            --report {output.report} \
            --min-bin {params.min_bin} &> {log}
        """


# Rule to analyze technical replicates
# Compares replicate samples to assess genotyping error rates
# Uses the raw sorted VCF from the project's assembly (which should contain replicate pairs)
rule qc_replicate:
    input:
        vcf="results/{project}/filtered_data/{project}.raw_sorted.vcf.gz"
    output:
        report="results/{project}/qc/{project}.replicate_analysis.txt"
    log:
        "logs/{project}/replicate_analysis.log"
    benchmark:
        "benchmarks/{project}/replicate_analysis.txt"
    params:
        suffix=lambda wildcards: config["projects"][wildcards.project]["parameters"]["replicate_qc"].get("suffix", "_repl")
    conda:
        "../envs/vcfpy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/vcf_analyze_replicates.py \
            --input {input.vcf} \
            --suffix {params.suffix} \
            --output {output.report} &> {log}
        """


# Rule to run vcftools Hardy-Weinberg equilibrium test
# Tests for deviations from HWE which may indicate genotyping errors or selection
rule vcftools_hw:
    input:
        vcf="results/{project}/filtered_data/{project}.raw_sorted.vcf.gz"
    output:
        hwe="results/{project}/qc/{project}.hardy.hwe"
    log:
        "logs/{project}/vcftools_hw.log"
    benchmark:
        "benchmarks/{project}/vcftools_hw.txt"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
            --hardy \
            --out results/{wildcards.project}/qc/{wildcards.project}.hardy &> {log}
        """


# Rule to summarize Hardy-Weinberg equilibrium test results
rule qc_hw:
    input:
        hwe="results/{project}/qc/{project}.hardy.hwe",
        vcf="results/{project}/filtered_data/{project}.raw_sorted.vcf.gz"
    output:
        report="results/{project}/qc/{project}.hardy_report.txt"
    log:
        "logs/{project}/qc_hw.log"
    benchmark:
        "benchmarks/{project}/qc_hw.txt"
    conda:
        "../envs/pandas.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/summarize_hardy.py \
            --input {input.hwe} \
            --vcf {input.vcf} \
            --output {output.report} &> {log}
        """
