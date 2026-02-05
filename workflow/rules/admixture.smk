# Rule for ADMIXTURE
rule admixture:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        Q = "results/{project}/admixture/{project}.biallelic_snps_thinned.{k}.Q",
        P = "results/{project}/admixture/{project}.biallelic_snps_thinned.{k}.P"
    log:
        "logs/{project}/admixture.{k}.log"
    benchmark:
        "benchmarks/{project}/admixture.K{k}.txt"
    params:
        output_dir = "results/{project}/admixture/",
        prefix = lambda wildcards: f"{wildcards.project}.biallelic_snps_thinned"
    conda:
        "../envs/admixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["admixture"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["admixture"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["admixture"]["runtime"]
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
        lambda wildcards: expand(
            "logs/{project}/admixture.{k}.log",
            project=wildcards.project,
            k=config["projects"][wildcards.project]["parameters"]["k_values"]
        )
    output:
        "results/{project}/admixture/{project}.admixture.chooseK_results.txt"
    shell:
        """
        echo "Choose the lowest cross-validation (CV) error:" > {output}
        grep -h 'CV' {input} >> {output} || true
        """
