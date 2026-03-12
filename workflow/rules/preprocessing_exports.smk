## Rules for exporting data from VCF to other formats
# Rule to convert filtered VCF to PLINK format
# renames chromosomes to "0" for downstream compatibility with "--allow-extra-chr 0"
rule vcf_to_plink:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        bed = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.bed",
        bim = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.bim",
        fam = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.fam"
    log:
        "logs/{project}/vcf_to_plink.log"
    benchmark:
        "benchmarks/{project}/vcf_to_plink.txt"
    params:
        output_prefix = lambda wildcards: f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_thinned"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --out {params.output_prefix} \
              --allow-extra-chr 0 \
              --double-id &> {log}
        """

# Rule to convert filtered VCF to STRUCTURE format
rule vcf_to_structure:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        str = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.str"
    log:
        "logs/{project}/vcf_to_structure.log"
    benchmark:
        "benchmarks/{project}/vcf_to_structure.txt"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/vcf_to_structure.py \
            --vcf {input.vcf} \
            --out {output.str} &> {log}
        """

## Rules used for generating VCFs with different data thresholds - for checking PCA robustness to missing data
# Rule to filter VCF for missing data threshold
rule filter_missing_vcf:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        vcf = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.vcf.gz"
    log:
        "logs/{project}/filter_missing_vcf_{miss}.log"
    benchmark:
        "benchmarks/{project}/filter_missing_vcf_{miss}.txt"
    params:
        miss = lambda wildcards: wildcards.miss
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools view -i 'F_MISSING<{params.miss}' {input.vcf} -Oz -o {output.vcf} &> {log}
        """

# Rule to convert to PLINK VCFs filtered for missing data threshold
# renames chromosomes to "0" for downstream compatibility with "--allow-extra-chr 0"
rule missing_vcf_to_plink:
    input:
        vcf = rules.filter_missing_vcf.output.vcf
    output:
        bed = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.bed",
        bim = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.bim",
        fam = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.fam"
    log:
        "logs/{project}/missing_vcf_to_plink_{miss}.log"
    benchmark:
        "benchmarks/{project}/missing_vcf_to_plink_{miss}.txt"
    params:
        miss = lambda wildcards: wildcards.miss,
        output_prefix = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --out {params.output_prefix} \
              --allow-extra-chr 0 \
              --double-id &> {log}
        """

