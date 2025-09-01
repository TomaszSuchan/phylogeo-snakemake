import os

# Rule to compress and index VCF if necessary
rule bgzip_vcf:
    input:
        # Snakemake selects the first existing file from this list
        vcf=lambda wc: (
            config["ipyrad_prefix"] + ".vcf.gz"
            if os.path.exists(config["ipyrad_prefix"] + ".vcf.gz")
            else config["ipyrad_prefix"] + ".vcf"
        )
    output:
        vcf="filtered_data/original.vcf.gz",
        index="filtered_data/original.vcf.gz.csi"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        time="10:00"
    shell:
        """
        # Detect extension
        input_vcf="{input.vcf}"
        ext="${{input_vcf##*.}}"

        if [[ "$ext" == "vcf.gz" ]]; then
            echo "Input is gzipped, sorting with bcftools..."
            bcftools index -f {input.vcf}
            bcftools sort {input.vcf} -Oz -o {output.vcf}
        else
            echo "Input is uncompressed, compressing and sorting..."
            bcftools sort {input.vcf} -Oz -o {output.vcf}
        fi

        # Index the output
        bcftools index -f {output.vcf}
        """

# Rule to select only biallelic SNPs with MAC>1 from VCF
rule select_biallelic_snps:
    input:
        vcf = rules.bgzip_vcf.output.vcf,
        index = rules.bgzip_vcf.output.index
    output:
        biallelic_vcf = "filtered_data/biallelic_snps.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb = 4000,
        time = "10:00"
    shell:
        """
        bcftools view -i 'MAC > 1' -m2 -M2 \
        -v snps {input.vcf} \
        -Oz -o {output.biallelic_vcf}
        """

# Rule to thin VCF using the custom script
rule thin_vcf:
    input:
        vcf = rules.select_biallelic_snps.output.biallelic_vcf
    output:
        thinned_vcf = "filtered_data/biallelic_snps_thinned.vcf.gz"
    params:
        min_coverage = config["vcf_thinning"].get("min_coverage", 0),
        method = config["vcf_thinning"].get("method", "max_coverage"),
        ties = config["vcf_thinning"].get("ties", "random"),
        ns_tag = config["vcf_thinning"].get("ns_tag", "NS"),
        id_pattern = config["vcf_thinning"].get("id_pattern", r"loc(\d+)_")
    conda:
        "../envs/vcfpy.yaml"
    threads: 1
    resources:
        mem_mb = 4000,
        time = "10:00"
    shell:
        """
        python workflow/scripts/thin_ipyrad_vcf.py \
            --vcf {input.vcf} \
            --out {output.thinned_vcf} \
            --min-coverage {params.min_coverage} \
            --method {params.method} \
            --ties {params.ties} \
            --ns-tag {params.ns_tag} \
            --id-pattern '{params.id_pattern}'
        """

# Rule to convert thinned VCF to PLINK format
rule vcf_to_plink:
    input:
        vcf = rules.thin_vcf.output.thinned_vcf
    output:
        bed = "filtered_data/biallelic_snps_thinned.bed",
        bim = "filtered_data/biallelic_snps_thinned.bim", 
        fam = "filtered_data/biallelic_snps_thinned.fam"
    params:
        output_prefix = "filtered_data/biallelic_snps_thinned"
    conda:
        "../envs/plink.yaml"
    threads: 1
    resources:
        mem_mb = 4000,
        time = "10:00"
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --out {params.output_prefix} \
              --allow-extra-chr \
              --double-id
        """

# Rule to convert thinned VCF to STRUCTURE format
rule vcf_to_structure:
    input:
        vcf = rules.thin_vcf.output.thinned_vcf
    output:
        str = "filtered_data/biallelic_snps_thinned.str"
    conda:
        "../envs/vcfpy.yaml"
    threads: 1
    resources:
        mem_mb = 4000,
        time = "10:00"
    shell:
        """
        python workflow/scripts/vcf_to_structure.py \
            --vcf {input.vcf} \
            --out {output.str}
        """

# Rule to filter VCF for missing data threshold and convert to PLINK
rule filter_missing_vcf:
    input:
        vcf = rules.thin_vcf.output.thinned_vcf
    output:
        vcf = "filtered_data/biallelic_snps_thinned_mincov{mincov}.vcf.gz"
    params:
        mincov = lambda wc: "{:.6f}".format(max(0.0, min(1.0, 1.0 - float(wc.mincov))))
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -i 'F_MISSING<{params.mincov}' {input.vcf} -Oz -o {output.vcf}
        """

# Rule to filter VCF for missing data threshold and convert to PLINK
rule missing_vcf_to_plink:
    input:
        vcf = rules.filter_missing_vcf.output.vcf
    output:
        bed = "filtered_data/biallelic_snps_thinned_mincov{mincov}.bed",
        bim = "filtered_data/biallelic_snps_thinned_mincov{mincov}.bim",
        fam = "filtered_data/biallelic_snps_thinned_mincov{mincov}.fam"
    params:
        mincov = lambda wildcards: wildcards.mincov,
        output_prefix = "filtered_data/biallelic_snps_thinned_mincov{mincov}"
    conda:
        "../envs/plink.yaml"
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --out {params.output_prefix} \
              --allow-extra-chr \
              --double-id
        """

# Rule to run complete preprocessing pipeline
rule run_preprocessing:
    input:
        rules.vcf_to_plink.output.bed,
        rules.vcf_to_structure.output.str,
        expand("filtered_data/biallelic_snps_thinned_mincov{mincov}.bed",
               mincov=config["PCA"]["mincov_thresholds"]),
        expand("filtered_data/biallelic_snps_thinned_mincov{mincov}.bim",
               mincov=config["PCA"]["mincov_thresholds"]),
        expand("filtered_data/biallelic_snps_thinned_mincov{mincov}.fam",
               mincov=config["PCA"]["mincov_thresholds"])
