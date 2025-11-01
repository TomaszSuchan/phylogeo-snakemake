# Rule for ADMIXTURE
rule admixture:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        Q = "{project}/admixture/biallelic_snps_thinned.{k}.Q",
        P = "{project}/admixture/biallelic_snps_thinned.{k}.P"
    log:
        "{project}/logs/admixture.{k}.log"
    benchmark:
        "{project}/benchmarks/admixture.K{k}.txt"
    params:
        output_dir = "{project}/admixture/",
        prefix = "biallelic_snps_thinned"
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
            "{project}/logs/admixture.{k}.log",
            project=wildcards.project,
            k=config["projects"][wildcards.project]["parameters"]["k_values"]
        )
    output:
        "{project}/admixture/chooseK_results.txt"
    shell:
        """
        echo "Choose the lowest cross-validation (CV) error:" > {output}
        grep -h 'CV' {input} >> {output} || true
        """
