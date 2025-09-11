# Rule for fastStructure
rule faststructure:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        meanQ = "results/faststructure/faststructure.{k}.meanQ",
        meanP = "results/faststructure/faststructure.{k}.meanP"
    params:
        input_prefix = rules.vcf_to_plink.params.output_prefix,
        output_prefix = "results/faststructure/faststructure",
        tol = config["faststructure"].get("tol", "10e-6"),
        prior = config["faststructure"].get("prior", "simple")
    conda:
        "../envs/faststructure.yaml"
    threads: config["resources"]["faststructure"]["threads"]
    resources:
        mem_mb = config["resources"]["faststructure"]["mem_mb"],
        time = config["resources"]["faststructure"]["runtime"]
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

# Rule to choose optimal K for fastStructure
rule faststructure_chooseK:
    input:
        expand("results/faststructure/faststructure.{k}.meanQ", k=config["k_values"])
    output:
        "results/faststructure/chooseK_results.txt"
    params:
        input_prefix =  "results/faststructure/faststructure",
    conda:
        "../envs/faststructure.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        time = config["resources"]["default"]["runtime"]
    shell:
        """
        chooseK.py \
            --input={params.input_prefix} \
            > {output}
        """

# Rule to run fastStructure for all K values
rule run_faststructure:
    input:
        rules.faststructure.output.meanQ,
        rules.faststructure.output.meanP,
        rules.faststructure_chooseK.output

