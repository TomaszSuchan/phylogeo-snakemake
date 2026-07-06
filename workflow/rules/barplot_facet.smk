# Facet panels combining per-K ancestry barplot RDS objects (K >= 2).

def _barplot_facet_plot_params(wildcards, method_label):
    """Layout for multi-K barplot facet panels (single column)."""
    mm = config["projects"][wildcards.project]["parameters"]["mapmixture"]
    facet = config["projects"][wildcards.project]["parameters"]["barplot_facet_plot"]
    return {
        "method_label": method_label,
        "label_width": facet["label_width"],
        "panel_gap": facet["panel_gap"],
        "legend_pad": facet.get("legend_pad", 0.4),
        "flip_axis": mm.get("flip_axis", False),
    }


rule plot_structure_barplot_facet:
    input:
        unpack(_structure_barplot_facet_inputs),
    output:
        pdf="results/{project}/structure/plots/{project}.structure.barplot-facet.pdf",
        rds="results/{project}/structure/plots/{project}.structure.barplot-facet.rds",
    params:
        lambda wildcards: _barplot_facet_plot_params(wildcards, "STRUCTURE"),
    log:
        "logs/{project}/plot_structure_barplot_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet.R"


rule plot_faststructure_barplot_facet:
    input:
        unpack(_faststructure_barplot_facet_inputs),
    output:
        pdf="results/{project}/faststructure/plots/{project}.faststructure.barplot-facet.pdf",
        rds="results/{project}/faststructure/plots/{project}.faststructure.barplot-facet.rds",
    params:
        lambda wildcards: _barplot_facet_plot_params(wildcards, "fastStructure"),
    log:
        "logs/{project}/plot_faststructure_barplot_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet.R"


rule plot_admixture_barplot_facet:
    input:
        unpack(_admixture_barplot_facet_inputs),
    output:
        pdf="results/{project}/admixture/plots/{project}.admixture.barplot-facet.pdf",
        rds="results/{project}/admixture/plots/{project}.admixture.barplot-facet.rds",
    params:
        lambda wildcards: _barplot_facet_plot_params(wildcards, "ADMIXTURE"),
    log:
        "logs/{project}/plot_admixture_barplot_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet.R"


rule plot_snmf_barplot_facet:
    input:
        unpack(_snmf_barplot_facet_inputs),
    output:
        pdf="results/{project}/snmf/plots/{project}.snmf.barplot-facet.pdf",
        rds="results/{project}/snmf/plots/{project}.snmf.barplot-facet.rds",
    params:
        lambda wildcards: _barplot_facet_plot_params(wildcards, "sNMF"),
    log:
        "logs/{project}/plot_snmf_barplot_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet.R"


rule plot_tess3_barplot_facet:
    input:
        unpack(_tess3_barplot_facet_inputs),
    output:
        pdf="results/{project}/tess3/plots/{project}.tess3.barplot-facet.pdf",
        rds="results/{project}/tess3/plots/{project}.tess3.barplot-facet.rds",
    params:
        lambda wildcards: _barplot_facet_plot_params(wildcards, "tess3"),
    log:
        "logs/{project}/plot_tess3_barplot_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet.R"


rule plot_alstructure_barplot_facet:
    input:
        unpack(_alstructure_barplot_facet_inputs),
    output:
        pdf="results/{project}/alstructure/plots/{project}.alstructure.barplot-facet.pdf",
        rds="results/{project}/alstructure/plots/{project}.alstructure.barplot-facet.rds",
    params:
        lambda wildcards: _barplot_facet_plot_params(wildcards, "ALStructure"),
    log:
        "logs/{project}/plot_alstructure_barplot_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet.R"


rule plot_construct_barplot_facet:
    input:
        unpack(_construct_barplot_facet_inputs),
    output:
        pdf="results/{project}/construct/plots/{project}.construct.barplot-facet.pdf",
        rds="results/{project}/construct/plots/{project}.construct.barplot-facet.rds",
    params:
        lambda wildcards: _barplot_facet_plot_params(wildcards, "conStruct"),
    log:
        "logs/{project}/plot_construct_barplot_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet.R"
