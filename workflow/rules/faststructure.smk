# Rule for fastStructure
rule faststructure:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        meanQ = "{project}/faststructure/faststructure.{k}.meanQ",
        meanP = "{project}/faststructure/faststructure.{k}.meanP"
    log:
        "{project}/logs/faststructure.{k}.log"
    benchmark:
        "{project}/benchmarks/faststructure.{k}.txt"
    params:
        input_prefix = rules.vcf_to_plink.params.output_prefix,
        output_prefix = "{project}/faststructure/faststructure",
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
            "{project}/faststructure/faststructure.{k}.meanQ", 
            project=wildcards.project,
            k=config["projects"][wildcards.project]["parameters"]["k_values"]
        )
    output:
        "{project}/faststructure/chooseK_results.txt"
    params:
        input_prefix =  "{project}/faststructure/faststructure",
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