import os

# Rule to sort input vcf
rule sort_vcf:
    input:
        # Use Snakemake's automatic file selection with multiple possible inputs
        vcf=lambda wildcards: (
            f"{config['ipyrad_prefix']}.vcf.gz"
            if os.path.exists(f"{config['ipyrad_prefix']}.vcf.gz")
            else f"{config['ipyrad_prefix']}.vcf"
        )
    output:
        vcf=config["analysis_name"] + "/filtered_data/raw_sorted.vcf.gz"
    log:
        config["analysis_name"] + "/logs/sort_vcf.log"
    conda:
        "../envs/vcftools.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
    shell:
        """
        vcf-sort {input.vcf} | bgzip -c > {output.vcf} &> {log}
        """

rule index_vcf:
    input:
        # Snakemake selects the first existing file from this list
        vcf=rules.sort_vcf.output.vcf
    output:
        index=config["analysis_name"] + "/filtered_data/raw_sorted.vcf.gz.csi"
    log:
        config["analysis_name"] + "/logs/index_vcf.log"
    conda:
        "../envs/bcftools.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
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
        biallelic_vcf = config["analysis_name"] + "/filtered_data/biallelic_snps.vcf.gz"
    log:
        config["analysis_name"] + "/logs/select_biallelic_snps.log
    conda:
        "../envs/bcftools.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
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
        thinned_vcf = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned.vcf.gz"
    log:
        config["analysis_name"] + "/logs/thin_vcf.log"
    params:
        min_coverage = config["vcf_thinning"].get("min_coverage", 0),
        method = config["vcf_thinning"].get("method", "max_coverage"),
        ties = config["vcf_thinning"].get("ties", "random"),
        ns_tag = config["vcf_thinning"].get("ns_tag", "NS"),
        id_pattern = config["vcf_thinning"].get("id_pattern", r"loc(\d+)_")
    conda:
        "../envs/vcfpy.yaml"
    threads: config["resources"]["default-long"]["threads"]
    resources:
        mem_mb = config["resources"]["default-long"]["mem_mb"],
        runtime = config["resources"]["default-long"]["runtime"]
    shell:
        """
        python workflow/scripts/thin_ipyrad_vcf.py \
            --vcf {input.vcf} \
            --out {output.thinned_vcf} \
            --min-coverage {params.min_coverage} \
            --method {params.method} \
            --ties {params.ties} \
            --ns-tag {params.ns_tag} \
            --id-pattern '{params.id_pattern}' &> {log}
        """

# Rule to convert thinned VCF to PLINK format
rule vcf_to_plink:
    input:
        vcf = rules.thin_vcf.output.thinned_vcf
    output:
        bed = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned.bed",
        bim = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned.bim",
        fam = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned.fam"
    log:
        config["analysis_name"] + "/logs/vcf_to_plink.log"
    params:
        output_prefix = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned"
    conda:
        "../envs/plink.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --out {params.output_prefix} \
              --allow-extra-chr \
              --double-id &> {log}
        """

# Rule to convert thinned VCF to STRUCTURE format
rule vcf_to_structure:
    input:
        vcf = rules.thin_vcf.output.thinned_vcf
    output:
        str = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned.str"
    log:
        config["analysis_name"] + "/logs/vcf_to_structure.log"
    conda:
        "../envs/vcfpy.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/vcf_to_structure.py \
            --vcf {input.vcf} \
            --out {output.str} &> {log}
        """

# Rule to filter VCF for missing data threshold and convert to PLINK
rule filter_missing_vcf:
    input:
        vcf = rules.thin_vcf.output.thinned_vcf
    output:
        vcf = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned_mincov{mincov}.vcf.gz"
    log:
        config["analysis_name"] + "/logs/filter_missing_vcf_{mincov}.log
    params:
        mincov = lambda wc: "{:.6f}".format(max(0.0, min(1.0, 1.0 - float(wc.mincov))))
    conda:
        "../envs/bcftools.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
    shell:
        """
        bcftools view -i 'F_MISSING<{params.mincov}' {input.vcf} -Oz -o {output.vcf} &> {log}
        """

# Rule to filter VCF for missing data threshold and convert to PLINK
rule missing_vcf_to_plink:
    input:
        vcf = rules.filter_missing_vcf.output.vcf
    output:
        bed = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned_mincov{mincov}.bed",
        bim = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned_mincov{mincov}.bim",
        fam = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned_mincov{mincov}.fam"
    log:
        config["analysis_name"] + "/logs/missing_vcf_to_plink_{mincov}.log"
    params:
        mincov = lambda wildcards: wildcards.mincov,
        output_prefix = "filtered_data/biallelic_snps_thinned_mincov{mincov}"
    conda:
        "../envs/plink.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --out {params.output_prefix} \
              --allow-extra-chr \
              --double-id &> {log}
        """
