"""
Rules for multivariate ordination analyses based on genetic distance matrices.
Includes PCoA on supported genetic distance matrices.
"""

def _get_pcoa_distance_input(wildcards):
    distance_inputs = {
        "euclidean": rules.euclidean_distance.output.dist,
        "kosman": rules.kosman_distance.output.dist,
        "avgsquared": rules.average_squared_genetic_difference.output.dist,
    }
    return distance_inputs[wildcards.distance]


rule pcoa_distance:
    """
    Perform PCoA (classical multidimensional scaling) on the requested
    distance matrix and output coordinates/eigenvalues in a PCA-like format.
    The output eigenvector file uses a first header row (ignored by the PCA
    plotting script) and has the first column with sample names so that
    it can be passed directly to the existing PCA plotting workflow.
    """
    input:
        dist=_get_pcoa_distance_input
    output:
        eigvecs="results/{project}/PCoA/{project}.PCoA_{distance}.eigvecs",
        eigvals="results/{project}/PCoA/{project}.PCoA_{distance}.eigvals"
    conda:
        "../envs/r-plot.yaml"
    log:
        "logs/{project}/pcoa_{distance}_distance.log",
    benchmark:
        "benchmarks/{project}/pcoa_{distance}_distance.txt",
    wildcard_constraints:
        distance="euclidean|kosman|avgsquared"
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
        eigvecs=rules.pcoa_distance.output.eigvecs,
        eigvals=rules.pcoa_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA_{distance}-PC{pc1}_PC{pc2}-{color_by}.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA_{distance}-PC{pc1}_PC{pc2}-{color_by}.rds"
    log:
        "logs/{project}/plot_pcoa_colored_{distance}_PC{pc1}_PC{pc2}_{color_by}.log"
    wildcard_constraints:
        distance="euclidean|kosman|avgsquared",
        color_by="(?!labeled|missing).*"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        color_by = lambda wildcards: wildcards.color_by,
        group_colors = lambda wildcards: _pca_plot_group_setting(
            wildcards.project, wildcards.color_by, "colors"
        ),
        plot_type = "colored",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
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
        eigvecs=rules.pcoa_distance.output.eigvecs,
        eigvals=rules.pcoa_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA_{distance}-PC{pc1}_PC{pc2}-labeled.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA_{distance}-PC{pc1}_PC{pc2}-labeled.rds"
    log:
        "logs/{project}/plot_pcoa_labeled_{distance}_PC{pc1}_PC{pc2}.log"
    wildcard_constraints:
        distance="euclidean|kosman|avgsquared"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "labeled",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
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
        eigvecs=rules.pcoa_distance.output.eigvecs,
        eigvals=rules.pcoa_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA_{distance}-PC{pc1}_PC{pc2}-missing.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA_{distance}-PC{pc1}_PC{pc2}-missing.rds"
    log:
        "logs/{project}/plot_pcoa_missing_{distance}_PC{pc1}_PC{pc2}.log"
    wildcard_constraints:
        distance="euclidean|kosman|avgsquared"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "missing",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
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
        eigvecs=rules.pcoa_distance.output.eigvecs,
        eigvals=rules.pcoa_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA_{distance}-facet-{color_by}.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA_{distance}-facet-{color_by}.rds"
    log:
        "logs/{project}/plot_pcoa_facet_colored_{distance}_{color_by}.log"
    wildcard_constraints:
        distance="euclidean|kosman|avgsquared",
        color_by="(?!labeled|missing).*"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        group_colors = lambda wildcards: _pca_plot_group_setting(
            wildcards.project, wildcards.color_by, "colors"
        ),
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "colored",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
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
        eigvecs=rules.pcoa_distance.output.eigvecs,
        eigvals=rules.pcoa_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA_{distance}-facet-labeled.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA_{distance}-facet-labeled.rds"
    log:
        "logs/{project}/plot_pcoa_facet_labeled_{distance}.log"
    wildcard_constraints:
        distance="euclidean|kosman|avgsquared"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "labeled",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
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
        eigvecs=rules.pcoa_distance.output.eigvecs,
        eigvals=rules.pcoa_distance.output.eigvals,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/PCoA/plots/{project}.PCoA_{distance}-facet-missing.pdf",
        rds="results/{project}/PCoA/plots/{project}.PCoA_{distance}-facet-missing.rds"
    log:
        "logs/{project}/plot_pcoa_facet_missing_{distance}.log"
    wildcard_constraints:
        distance="euclidean|kosman|avgsquared"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "missing",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    group: "plot_pcoa"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"




