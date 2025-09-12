# Rule for ADMIXTURE
rule admixture:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        Q = config["analysis_name"] + "/admixture/admixture.K{k}.Q",
        P = config["analysis_name"] + "/admixture/admixture.K{k}.P"
    log:
        config["analysis_name"] + "/logs/admixture.K{k}.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/admixture.K{k}.txt"
    params:
        output_dir = config["analysis_name"] + "/admixture/"
    conda:
        "../envs/admixture.yaml"
    threads: config["resources"]["admixture"]["threads"]
    resources:
        mem_mb = config["resources"]["admixture"]["mem_mb"],
        time = config["resources"]["admixture"]["runtime"]
    shell:
        """
        cd {params.output_dir}
        admixture --cv -j{threads} ../../{input.bed} {wildcards.k} &> {log}
        """