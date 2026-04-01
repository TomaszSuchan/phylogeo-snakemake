"""
Rules for phylogenetic network inference and visualization using NeighborNet.
NeighborNet is inferred from the Euclidean genetic distance matrix.
"""

rule neighbornet_euclidean:
    """
    Build a NeighborNet network from the Euclidean genetic distance matrix.
    """
    input:
        dist=rules.euclidean_distance.output.dist
    output:
        net="results/{project}/neighbornet/{project}.euclidean.neighbornet.rds"
    log:
        "logs/{project}/neighbornet_euclidean.log"
    benchmark:
        "benchmarks/{project}/neighbornet_euclidean.txt"
    conda:
        "../envs/neighbornet.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/build_neighbornet.R"


rule plot_neighbornet:
    """
    Plot NeighborNet with tip colors from a configured indpopdata column.
    """
    input:
        net=rules.neighbornet_euclidean.output.net,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        pdf="results/{project}/neighbornet/plots/{project}.euclidean.neighbornet-{color_by}.pdf",
        rds="results/{project}/neighbornet/plots/{project}.euclidean.neighbornet-{color_by}.rds"
    log:
        "logs/{project}/plot_neighbornet_{color_by}.log"
    params:
        color_by=lambda wildcards: wildcards.color_by,
        neighbornet_colors=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get(
            "neigbournet_colors",
            config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("colors", None)
        ),
        width=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("width", 12),
        height=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("height", 12),
        dpi=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("neighbornet", {}).get("dpi", 300)
    conda:
        "../envs/neighbornet.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_neighbornet.R"
