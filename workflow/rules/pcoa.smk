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
        eigvecs="results/{project}/PCoA/{project}.euclidean_pcoa.eigvecs",
        eigvals="results/{project}/PCoA/{project}.euclidean_pcoa.eigvals"
    conda:
        "../envs/r-plot.yaml"
    log:
        "logs/{project}/pcoa_euclidean_distance.log",
    benchmark:
        "benchmarks/{project}/pcoa_euclidean_distance.txt",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/pcoa_from_euclidean_distance.R"

rule plot_pcoa_colored:
    """
    Plot PCoA axes colored by population/metadata, reusing the PCA plotting script.
    """
    input:
        eigvecs=rules.pcoa_euclidean_distance.output.eigvecs,
        eigvals=rules.pcoa_euclidean_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA-PC{pc1}_PC{pc2}-{color_by}.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA-PC{pc1}_PC{pc2}-{color_by}.rds"
    log:
        "logs/{project}/plot_pcoa_colored_PC{pc1}_PC{pc2}_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        color_by = lambda wildcards: wildcards.color_by,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        plot_type = "colored"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    group: "plot_pcoa"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

rule plot_pcoa_labeled:
    """
    Plot PCoA axes with sample labels only (no coloring).
    """
    input:
        eigvecs=rules.pcoa_euclidean_distance.output.eigvecs,
        eigvals=rules.pcoa_euclidean_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA-PC{pc1}_PC{pc2}-labeled.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA-PC{pc1}_PC{pc2}-labeled.rds"
    log:
        "logs/{project}/plot_pcoa_labeled_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "labeled"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    group: "plot_pcoa"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

rule plot_pcoa_missing:
    """
    Plot PCoA axes colored by missing data per individual.
    """
    input:
        eigvecs=rules.pcoa_euclidean_distance.output.eigvecs,
        eigvals=rules.pcoa_euclidean_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA-PC{pc1}_PC{pc2}-missing.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA-PC{pc1}_PC{pc2}-missing.rds"
    log:
        "logs/{project}/plot_pcoa_missing_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "missing"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    group: "plot_pcoa"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

rule plot_pcoa_facet_colored:
    """
    Facet plot of all PCoA axis combinations, colored by metadata.
    """
    input:
        eigvecs=rules.pcoa_euclidean_distance.output.eigvecs,
        eigvals=rules.pcoa_euclidean_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA-facet-{color_by}.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA-facet-{color_by}.rds"
    log:
        "logs/{project}/plot_pcoa_facet_colored_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "colored"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    group: "plot_pcoa"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

rule plot_pcoa_facet_labeled:
    """
    Facet plot of all PCoA axis combinations with labels only.
    """
    input:
        eigvecs=rules.pcoa_euclidean_distance.output.eigvecs,
        eigvals=rules.pcoa_euclidean_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA-facet-labeled.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA-facet-labeled.rds"
    log:
        "logs/{project}/plot_pcoa_facet_labeled.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "labeled"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    group: "plot_pcoa"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

rule plot_pcoa_facet_missing:
    """
    Facet plot of all PCoA axis combinations colored by missing data.
    """
    input:
        eigvecs=rules.pcoa_euclidean_distance.output.eigvecs,
        eigvals=rules.pcoa_euclidean_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA-facet-missing.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA-facet-missing.rds"
    log:
        "logs/{project}/plot_pcoa_facet_missing.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "missing"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    group: "plot_pcoa"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"




