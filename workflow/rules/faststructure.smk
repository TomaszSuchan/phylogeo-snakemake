# Rule for fastStructure
rule faststructure:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        meanQ = "results/{project}/faststructure/{project}.faststructure.{k}.meanQ",
        meanP = "results/{project}/faststructure/{project}.faststructure.{k}.meanP"
    log:
        "logs/{project}/faststructure.{k}.log"
    benchmark:
        "benchmarks/{project}/faststructure.{k}.txt"
    params:
        input_prefix = lambda wildcards: f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_thinned",
        output_prefix = "results/{project}/faststructure/{project}.faststructure",
        tol = lambda wildcards:  config["projects"][wildcards.project]["parameters"]["faststructure"].get("tol", "10e-6"),
        prior = lambda wildcards: config["projects"][wildcards.project]["parameters"]["faststructure"].get("prior", "simple")
    conda:
        "../envs/faststructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["faststructure"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["faststructure"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["faststructure"]["runtime"]
    shell:
        """
        structure.py \
            -K {wildcards.k} \
            --input={params.input_prefix} \
            --output={params.output_prefix} \
            --tol={params.tol} \
            --prior={params.prior} \
            --cv=0 &> {log}
        """

# Rule to choose optimal K for fastStructure
rule faststructure_chooseK:
    input:
        lambda wildcards: expand(
            "results/{project}/faststructure/{project}.faststructure.{k}.meanQ",
            project=wildcards.project,
            k=config["projects"][wildcards.project]["parameters"]["k_values"]
        )
    output:
        "results/{project}/faststructure/{project}.faststructure.chooseK_results.txt"
    params:
        input_prefix =  "results/{project}/faststructure/{project}.faststructure",
    conda:
        "../envs/faststructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        chooseK.py \
            --input={params.input_prefix} \
            > {output}
        """