# Rule for structure
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

rule structure:
    input:
        ustr = config["ipyrad_prefix"] + ".ustr"
    output:
        meanQ = "faststructure/input_K{k}.meanQ",
        meanP = "faststructure/input_K{k}.meanP"
    params:
        input_prefix = rules.vcf_to_plink.params.output_prefix,
        output_prefix = "faststructure/input_K{k}",
        tol = config["faststructure"].get("faststructure_tol", "10e-6"),
        prior = config["faststructure"].get("faststructure_prior", "simple")
    conda:
        "../envs/faststructure.yaml"
    threads: 4
    resources:
        mem_mb = 8000,
        time = "2:00:00"
    shell:
        """
        structure.py \
            -K {wildcards.k} \
            --input={params.input_prefix} \
            --output={params.output_prefix} \
            --tol={params.tol} \
            --prior={params.prior} \
            --cv=0
        """