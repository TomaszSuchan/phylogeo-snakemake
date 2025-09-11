# Rule for ADMIXTURE
rule admixture:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        Q = "results/admixture/admixture.K{k}.Q",
        P = "results/admixture/admixture.K{k}.P"
    params:
        output_dir = "results/admixture/"
    conda:
        "../envs/admixture.yaml"
    threads: config["resources"]["admixture"]["threads"]
    resources:
        mem_mb = config["resources"]["admixture"]["mem_mb"],
        time = config["resources"]["admixture"]["runtime"]
    shell:
        """
        cd {params.output_dir}
        admixture --cv -j{threads} ../../{input.bed} {wildcards.k} | tee log{wildcards.k}.out
        """