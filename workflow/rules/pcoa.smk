"""
Rules for multivariate ordination analyses based on genetic distance matrices.
Currently includes PCoA on the Euclidean genetic distance matrix.
"""

rule pcoa_euclidean_distance:
    """
    Perform PCoA (classical multidimensional scaling) on the Euclidean
    distance matrix and output coordinates/eigenvalues in a PCA-like format.
    The output eigenvector file uses a first header row (ignored by the PCA
    plotting script) and has the first column with sample names so that
    it can be passed directly to the existing PCA plotting workflow.
    """
    input:
        dist=rules.euclidean_distance.output.dist
    output:
        eigvecs="results/{project}/gen_dist/{project}.euclidean_pcoa.eigvecs",
        eigvals="results/{project}/gen_dist/{project}.euclidean_pcoa.eigvals"
    conda:
        "../envs/r-plot.yaml"
    log:
        "logs/{project}/euclidean_pcoa.log",
    benchmark:
        "benchmarks/{project}/euclidean_pcoa.txt",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/pcoa_from_euclidean_distance.R"

