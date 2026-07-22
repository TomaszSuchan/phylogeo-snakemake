"""
Rules for phylogenetic network inference and visualization using NeighborNet.
NeighborNet is inferred from the p-distance genetic distance matrix via fastnntr.
"""

rule install_fastnntr:
    """
    Install fastnntr once into the Snakemake conda environment from GitHub.
    fastnntr is an extendr/Rust package (not on conda); the env provides rustc/cargo
    and this rule installs a pinned release before NeighborNet inference.
    """
    output:
        touch(".snakemake/fastnntr_installed")
    conda:
        "../envs/neighbornet.yaml"
    threads: config["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = config["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = config["parameters"]["resources"]["default-long"]["runtime"]
    log:
        "logs/install_fastnntr.log"
    shell:
        """
        set -euo pipefail
        Rscript --vanilla -e 'lib <- .libPaths()[1]; unlink(list.files(lib, pattern="^00LOCK-", full.names=TRUE), recursive=TRUE, force=TRUE); if (!requireNamespace("fastnntr", quietly=TRUE, lib.loc=lib)) remotes::install_github("rhysnewell/fast-nnt", subdir="fastnntr", ref="v0.2.5", upgrade="never", dependencies=TRUE, build_vignettes=FALSE, lib=lib); if (!requireNamespace("fastnntr", quietly=TRUE, lib.loc=lib)) stop("Failed to install required R package: fastnntr")' &> {log}
        touch {output}
        """

rule neighbornet_pdistance:
    """
    Build a NeighborNet network from the p-distance genetic distance matrix.
    """
    input:
        dist=rules.p_distance.output.dist,
        install=rules.install_fastnntr.output
    output:
        net="results/{project}/neighbornet/{project}.pdistance.neighbornet.rds"
    log:
        "logs/{project}/neighbornet_pdistance.log"
    benchmark:
        "benchmarks/{project}/neighbornet_pdistance.txt"
    params:
        ordering_method=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("ordering_method", "closest-pair"),
        inference_method=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("inference_method", "active-set"),
    conda:
        "../envs/neighbornet.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["neighbornet"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["neighbornet"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["neighbornet"]["runtime"]
    script:
        "../scripts/build_neighbornet.R"


rule plot_neighbornet:
    """
    Plot NeighborNet with tip colors from a configured indpopdata column.
    """
    input:
        net=rules.neighbornet_pdistance.output.net,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        pdf_with_tip_labels="results/{project}/neighbornet/plots/{project}.pdistance.neighbornet-{color_by}.with_tip_labels.pdf",
        pdf_no_tip_labels="results/{project}/neighbornet/plots/{project}.pdistance.neighbornet-{color_by}.no_tip_labels.pdf",
        pdf_phangorn="results/{project}/neighbornet/plots/{project}.pdistance.neighbornet-{color_by}.phangorn.pdf",
        rds="results/{project}/neighbornet/plots/{project}.pdistance.neighbornet-{color_by}.rds",
    log:
        "logs/{project}/plot_neighbornet_{color_by}.log"
    params:
        color_by=lambda wildcards: wildcards.color_by,
        group_colors=lambda wildcards: _neighbornet_group_setting(
            wildcards.project, wildcards.color_by, "colors"
        ),
        width=lambda wildcards: _fig_cm_to_in(config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("width"), 30.48),
        height=lambda wildcards: _fig_cm_to_in(config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("height"), 30.48),
        dpi=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("dpi", 300),
        linewidth=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("linewidth", 0.5),
    conda:
        "../envs/neighbornet.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_neighbornet.R"
