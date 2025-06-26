# Rule for fastStructure
rule faststructure:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        meanQ = "faststructure/input_K{k}.meanQ",
        meanP = "faststructure/input_K{k}.meanP"
    params:
        input_prefix = rules.vcf_to_plink.params.output_prefix,
        output_prefix = "faststructure/input_K{k}",
        tol = config["faststructure"].get("faststructure_tol", "10e-6"),
        prior = config["faststructure"].get("faststructure_prior", "simple")
    conda:
        "../envs/faststructure.yaml"
    threads: 4
    resources:
        mem_mb = 8000,
        time = "2:00:00"
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
        expand("faststructure/input_K{k}.meanQ", k=config["k_values"])
    output:
        "faststructure/chooseK_results.txt"
    params:
        input_prefix = "faststructure/input"
    conda:
        "../envs/faststructure.yaml"
    threads: 1
    resources:
        mem_mb = 2000,
        time = "10:00"
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

