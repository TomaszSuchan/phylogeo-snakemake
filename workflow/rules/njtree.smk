"""
Neighbour-joining tree with rapidNJ from the filtered .loci alignment.
"""

rule njtree:
    input:
        fasta=rules.loci_to_alignment.output.fasta,
    output:
        treefile="results/{project}/njtree/{project}.njtree.treefile",
    log:
        "logs/{project}/njtree.log",
    benchmark:
        "benchmarks/{project}/njtree.txt",
    params:
        distance_model=lambda wildcards: config["projects"][wildcards.project]["parameters"]["njtree"].get("distance_model", "kim"),
        bootstrap_flag=lambda wildcards: (
            f"-b {int(config['projects'][wildcards.project]['parameters']['njtree'].get('bootstraps', 1000))}"
            if int(config["projects"][wildcards.project]["parameters"]["njtree"].get("bootstraps", 1000)) > 0
            else ""
        ),
    conda:
        "../envs/rapidnj.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["njtree"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["njtree"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["njtree"]["runtime"],
    shell:
        """
        # -a is the sequence distance model (jc/kim). Do not use -d here: that is HDD cache dir.
        rapidnj {input.fasta} -i fa -a {params.distance_model} -c {threads} \
            -x {output.treefile} {params.bootstrap_flag} > {log} 2>&1
        """


rule plot_njtree:
    input:
        treefile=rules.njtree.output.treefile,
    output:
        unrooted_pdf="results/{project}/njtree/plots/{project}.tree_plot.unrooted.pdf",
        unrooted_rds="results/{project}/njtree/plots/{project}.tree_plot.unrooted.rds",
        unrooted_tree="results/{project}/njtree/{project}.njtree.unrooted.treefile",
    benchmark:
        "benchmarks/{project}/plot_njtree.txt",
    params:
        support_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"]["njtree"].get("support_threshold", 70),
    conda:
        "../envs/r-ggtree.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("mem_mb", 8000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("runtime", 10),
    script:
        "../scripts/plot_nj_tree.R"
