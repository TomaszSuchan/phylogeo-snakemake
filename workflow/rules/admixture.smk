# Rule for ADMIXTURE
rule admixture:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        Q = config["analysis_name"] + "/admixture/biallelic_snps_thinned.{k}.Q",
        P = config["analysis_name"] + "/admixture/biallelic_snps_thinned.{k}.P"
    log:
        config["analysis_name"] + "/logs/admixture.{k}.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/admixture.K{k}.txt"
    params:
        output_dir = config["analysis_name"] + "/admixture/",
        prefix = lambda wildcards, input: input.bed.rsplit("/",1)[-1].replace(".bed","")
    conda:
        "../envs/admixture.yaml"
    threads: config["resources"]["admixture"]["threads"]
    resources:
        mem_mb = config["resources"]["admixture"]["mem_mb"],
        runtime = config["resources"]["admixture"]["runtime"]
    shell:
        """
        mkdir -p {params.output_dir}
        admixture --cv -j{threads} {input.bed} {wildcards.k} 2>&1 | tee {log}
        mv {params.prefix}.{wildcards.k}.Q {params.output_dir}/
        mv {params.prefix}.{wildcards.k}.P {params.output_dir}/
        """

# Rule to choose optimal K for admixture
rule admixture_chooseK:
    input:
        expand(config["analysis_name"] + "/logs/admixture.{k}.log", k=config["k_values"])
    output:
        config["analysis_name"] + "/admixture/chooseK_results.txt"
    shell:
        """
        echo "Choose the lowest cross-validation (CV) error:" > {output}
        grep -h 'CV' {input} >> {output} || true
        """