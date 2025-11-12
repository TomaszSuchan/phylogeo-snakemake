rule iqtree:
    input:
        phy=lambda wildcards: f"{config['projects'][wildcards.project]['ipyrad_prefix']}.phy"
    output:
        treefile="results/{project}/iqtree/{project}.iqtree.treefile",
        iqtree="results/{project}/iqtree/{project}.iqtree.iqtree",
        log="results/{project}/iqtree/{project}.iqtree.log"
    benchmark:
        "benchmarks/{project}/iqtree.txt"
    params:
        prefix="results/{project}/iqtree/{project}.iqtree",
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
        phy=lambda wildcards: f"{config['projects'][wildcards.project]['ipyrad_prefix']}.phy"
    output:
        treefile="results/{project}/iqtree_robust/{project}.iqtree_robust.treefile",
        iqtree="results/{project}/iqtree_robust/{project}.iqtree_robust.iqtree",
        log="results/{project}/iqtree_robust/{project}.iqtree_robust.log"
    benchmark:
        "benchmarks/{project}/iqtree_robust.txt"
    params:
        prefix="results/{project}/iqtree_robust/{project}.iqtree_robust",
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
        pdf="results/{project}/iqtree/{project}.tree_plot.pdf",
        rds="results/{project}/iqtree/{project}.tree_plot.rds"
    benchmark:
        "benchmarks/{project}/plot_tree.txt"
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
        treefile="results/{project}/iqtree_robust/{project}.iqtree_robust.treefile"
    output:
        pdf="results/{project}/iqtree_robust/{project}.tree_plot.pdf",
        rds="results/{project}/iqtree_robust/{project}.tree_plot.rds"
    benchmark:
        "benchmarks/{project}/plot_tree_robust.txt"
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
