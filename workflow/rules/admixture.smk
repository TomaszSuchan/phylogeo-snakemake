# Rule for ADMIXTURE
rule admixture:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        Q = "admixture/input.{k}.Q",
        P = "admixture/input.{k}.P"
    params:
        output_dir = "admixture/"
    conda:
        "../envs/admixture.yaml"
    threads: 4
    resources:
        mem_mb = 8000,
        time = "1:00:00"
    shell:
        """
        cd {params.output_dir}
        admixture -j{threads} ../{input.bed} {wildcards.k} | tee log{wildcards.k}.out
        """

# Rule to run ADMIXTURE for all K values  
rule run_admixture:
    input:
        rules.admixture.output.Q,
        rules.admixture.output.P
