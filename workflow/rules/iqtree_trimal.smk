rule N_to_gaps:
    input:
        phy=lambda wildcards: f"{config['projects'][wildcards.project]['ipyrad_prefix']}.phy"
    output:
        gapped_phy="results/{project}/filtered_data/{project}.gapped.phy"
    benchmark:
        "benchmarks/{project}/N_to_gaps.txt"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("mem_mb", 4000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("runtime", 10)
    shell:
        """
        sed 's/N/-/g' {input.phy} > {output.gapped_phy}
        """

rule trimal:
    input:
        phy=rules.N_to_gaps.output.gapped_phy
    output:
        trimmed_phy="results/{project}/filtered_data/{project}.trimmed.phy"
    benchmark:
        "benchmarks/{project}/trimal.txt"
    params:
        gap_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"]["trimal"].get("gt", 0.75)
    conda:
        "../envs/trimal.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("mem_mb", 4000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("runtime", 10)
    shell:
        """
        trimal -in {input.phy} -out {output.trimmed_phy} -gt {params.gap_threshold}
        """

rule iqtree:
    input:
        phy=rules.trimal.output.trimmed_phy
    output:
        treefile="results/{project}/iqtree_trimal/{project}.iqtree.treefile",
        iqtree="results/{project}/iqtree_trimal/{project}.iqtree.iqtree",
        log="results/{project}/iqtree_trimal/{project}.iqtree.log"
    benchmark:
        "benchmarks/{project}/iqtree_trimal.txt"
    params:
        prefix="results/{project}/iqtree_trimal/{project}.iqtree",
        model=lambda wildcards: config["projects"][wildcards.project]["parameters"]["iqtree"].get("model", "MFP"),
        bootstraps=lambda wildcards: config["projects"][wildcards.project]["parameters"]["iqtree"].get("bootstraps", "1000"),
        outgroup=lambda wildcards: (f"-o {config['projects'][wildcards.project]['parameters']['iqtree']['outgroup']}"
                                    if "outgroup" in config["projects"][wildcards.project]["parameters"].get("iqtree", {})
                                    else "")
    conda:
        "../envs/iqtree.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["iqtree"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["iqtree"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["iqtree"]["runtime"]
    shell:
        """
        iqtree -s {input.phy} -m {params.model} -B {params.bootstraps} \
               {params.outgroup} -T {threads} -pre {params.prefix}
        """

rule iqtree_robust:
    input:
        phy=rules.trimal.output.trimmed_phy
    output:
        treefile="results/{project}/iqtree_trimal_robust/{project}.iqtree_robust.treefile",
        iqtree="results/{project}/iqtree_trimal_robust/{project}.iqtree_robust.iqtree",
        log="results/{project}/iqtree_trimal_robust/{project}.iqtree_robust.log"
    benchmark:
        "benchmarks/{project}/iqtree_trimal_robust.txt"
    params:
        prefix="results/{project}/iqtree_trimal_robust/{project}.iqtree_robust",
        model=lambda wildcards: config["projects"][wildcards.project]["parameters"]["iqtree"].get("model", "MFP"),
        robust_phy=lambda wildcards: config['projects'][wildcards.project]["parameters"]["iqtree"].get("robust-phy", "0.95"),
        bootstraps=lambda wildcards: config["projects"][wildcards.project]["parameters"]["iqtree"].get("bootstraps", "1000"),
        outgroup=lambda wildcards: (f"-o {config['projects'][wildcards.project]['parameters']['iqtree']['outgroup']}"
                                    if "outgroup" in config["projects"][wildcards.project]["parameters"].get("iqtree", {})
                                    else "")
    conda:
        "../envs/iqtree.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["iqtree"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["iqtree"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["iqtree"]["runtime"]
    shell:
        """
        iqtree -s {input.phy} -m {params.model} -B {params.bootstraps} \
               {params.outgroup} -T {threads} -pre {params.prefix} --robust-phy {params.robust_phy}
        """

rule plot_tree:
    input:
        treefile=rules.iqtree.output.treefile
    output:
        pdf="results/{project}/iqtree_trimal/{project}.tree_plot.pdf",
        rds="results/{project}/iqtree_trimal/{project}.tree_plot.rds"
    benchmark:
        "benchmarks/{project}/iqtree_trimal_plot.txt"
    params:
        support_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"]["iqtree"].get("support_threshold", 70)
    conda:
        "../envs/r-ggtree.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("mem_mb", 8000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("runtime", 10)
    script:
        "../scripts/plot_tree.R"


rule plot_tree_robust:
    input:
        treefile=rules.iqtree_robust.output.treefile
    output:
        pdf="results/{project}/iqtree_trimal_robust/{project}.tree_plot.pdf",
        rds="results/{project}/iqtree_trimal_robust/{project}.tree_plot.rds"
    benchmark:
        "benchmarks/{project}/iqtree_trimal_robust_plot.txt"
    params:
        support_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"]["iqtree"].get("support_threshold", 70)
    conda:
        "../envs/r-ggtree.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("mem_mb", 8000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("runtime", 10)
    script:
        "../scripts/plot_tree.R"
