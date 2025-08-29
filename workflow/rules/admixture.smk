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
    threads: 4
    resources:
        mem_mb = 8000,
        time = "1:00:00"
    shell:
        """
        cd {params.output_dir}
        admixture --cv -j{threads} ../../{input.bed} {wildcards.k} | tee log{wildcards.k}.out
        """

# Rule to run ADMIXTURE for all K values  
rule run_admixture:
    input:
        expand("results/admixture/admixture.K{k}.Q", k=config["k_values"]),
        expand("results/admixture/admixture.K{k}.P", k=config["k_values"])
