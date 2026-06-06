"""
Rules for phylogenetic network inference and visualization using NeighborNet.
NeighborNet is inferred from the p-distance genetic distance matrix.
"""

rule neighbornet_pdistance:
    """
    Build a NeighborNet network from the p-distance genetic distance matrix.
    """
    input:
        dist=rules.p_distance.output.dist
    output:
        net="results/{project}/neighbornet/{project}.pdistance.neighbornet.rds"
    log:
        "logs/{project}/neighbornet_pdistance.log"
    benchmark:
        "benchmarks/{project}/neighbornet_pdistance.txt"
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
        neighbornet_colors=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get(
            "neigbournet_colors", None
        ),
        width=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("width", 12),
        height=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("height", 12),
        dpi=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("dpi", 300),
    conda:
        "../envs/neighbornet.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_neighbornet.R"
