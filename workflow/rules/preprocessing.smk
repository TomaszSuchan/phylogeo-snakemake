import os

# Rule to sort input vcf
rule sort_vcf:
    input:
      vcf=lambda wildcards: (
          f"{config['projects'][wildcards.project]['ipyrad_prefix']}.vcf.gz"
          if os.path.exists(f"{config['projects'][wildcards.project]['ipyrad_prefix']}.vcf.gz")
          else f"{config['projects'][wildcards.project]['ipyrad_prefix']}.vcf"
      )
    output:
        vcf="results/{project}/filtered_data/{project}.raw_sorted.vcf.gz"
    log:
        "logs/{project}/sort_vcf.log"
    benchmark:
        "benchmarks/{project}/sort_vcf.txt"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcf-sort {input.vcf} 2> {log} | bgzip -c > {output.vcf}
        """

rule index_vcf:
    input:
        # Snakemake selects the first existing file from this list
        vcf=rules.sort_vcf.output.vcf
    output:
        index="results/{project}/filtered_data/{project}.raw_sorted.vcf.gz.csi"
    log:
        "logs/{project}/index_vcf.log"
    benchmark:
        "benchmarks/{project}/index_vcf.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools index -f {input.vcf} &> {log}
        """

# Rule to select only biallelic SNPs with MAC>1 from VCF
rule select_biallelic_snps:
    input:
        vcf = rules.sort_vcf.output.vcf,
        index = rules.index_vcf.output.index
    output:
        biallelic_vcf = "results/{project}/filtered_data/{project}.biallelic_snps.vcf.gz"
    log:
        "logs/{project}/select_biallelic_snps.log"
    benchmark:
        "benchmarks/{project}/select_biallelic_snps.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools view -i 'MAC > 1' -m2 -M2 \
        -v snps {input.vcf} \
        -Oz -o {output.biallelic_vcf} &> {log}
        """

# Rule to thin VCF using the custom script
rule thin_vcf:
    input:
        vcf = rules.select_biallelic_snps.output.biallelic_vcf
    output:
        vcf = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.vcf.gz"
    log:
        "logs/{project}/thin_vcf.log"
    benchmark:
        "benchmarks/{project}/thin_vcf.txt"
    params:
        min_coverage = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("min_coverage", 0),
        method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("method", "max_coverage"),
        ties = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("ties", "random"),
        ns_tag = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("ns_tag", "NS"),
        id_pattern = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("id_pattern", r"loc(\d+)_")
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/thin_ipyrad_vcf.py \
            --vcf {input.vcf} \
            --out {output.vcf} \
            --min-coverage {params.min_coverage} \
            --method {params.method} \
            --ties {params.ties} \
            --ns-tag {params.ns_tag} \
            --id-pattern '{params.id_pattern}' &> {log}
        """

# Rule to convert thinned VCF to PLINK format
# renames chromosomes to "0" for downstream compatibility with "--allow-extra-chr 0"
rule vcf_to_plink:
    input:
        vcf = rules.thin_vcf.output.vcf
    output:
        bed = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.bed",
        bim = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.bim",
        fam = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.fam"
    log:
        "logs/{project}/vcf_to_plink.log"
    benchmark:
        "benchmarks/{project}/vcf_to_plink.txt"
    params:
        output_prefix = "results/{project}/filtered_data/{project}.biallelic_snps_thinned"
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

# Rule to convert thinned VCF to STRUCTURE format
rule vcf_to_structure:
    input:
        vcf = rules.thin_vcf.output.vcf
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

# Rule to filter VCF for missing data threshold
rule filter_missing_vcf:
    input:
        vcf = rules.thin_vcf.output.vcf
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
