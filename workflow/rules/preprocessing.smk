# Rule to select only biallelic SNPs with MAC>1 from VCF
rule select_biallelic_snps:
    input:
        vcf = config["ipyrad_prefix"] + ".vcf"
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
        bcftools view -i 'MAC <= 1' -m2 -M2 \
        -v snps {input.vcf} \
        -Oz -o {output.biallelic_vcf}
        """

# Rule to thin VCF using ipyrad script
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
        plink2 --vcf {input.vcf} \
               --make-bed \
               --out {params.output_prefix} \
               --allow-extra-chr \
               --double-id
        """

# Rule to run complete preprocessing pipeline
rule run_preprocessing:
    input:
        rules.vcf_to_plink.output
