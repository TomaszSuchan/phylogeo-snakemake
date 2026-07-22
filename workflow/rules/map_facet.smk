# Facet panels combining per-K ancestry map panel RDS objects (K >= 2).
# patchwork::wrap_plots on pre-rendered grobs; choose-K plots are separate outputs.

def _map_facet_plot_params(wildcards, method_label):
    """Layout for multi-K map facet panels (page size follows map_background)."""
    bg = config["projects"][wildcards.project]["parameters"]["map_background"]
    facet = config["projects"][wildcards.project]["parameters"].get("map_facet_plot", {})
    return {
        "method_label": method_label,
        "ncol": facet.get("ncol", 2),
        "panel_width": _fig_cm_to_in(bg.get("width"), 25.4),
        "panel_height": _fig_cm_to_in(bg.get("height"), 20.32),
    }


rule plot_structure_map_facet:
    input:
        unpack(_structure_map_facet_inputs),
    output:
        pdf="results/{project}/structure/plots/{project}.structure.map-facet.pdf",
        rds="results/{project}/structure/plots/{project}.structure.map-facet.rds",
    params:
        lambda wildcards: _map_facet_plot_params(wildcards, "STRUCTURE"),
    log:
        "logs/{project}/plot_structure_map_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_map_facet.R"


rule plot_faststructure_map_facet:
    input:
        unpack(_faststructure_map_facet_inputs),
    output:
        pdf="results/{project}/faststructure/plots/{project}.faststructure.map-facet.pdf",
        rds="results/{project}/faststructure/plots/{project}.faststructure.map-facet.rds",
    params:
        lambda wildcards: _map_facet_plot_params(wildcards, "fastStructure"),
    log:
        "logs/{project}/plot_faststructure_map_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_map_facet.R"


rule plot_admixture_map_facet:
    input:
        unpack(_admixture_map_facet_inputs),
    output:
        pdf="results/{project}/admixture/plots/{project}.admixture.map-facet.pdf",
        rds="results/{project}/admixture/plots/{project}.admixture.map-facet.rds",
    params:
        lambda wildcards: _map_facet_plot_params(wildcards, "ADMIXTURE"),
    log:
        "logs/{project}/plot_admixture_map_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_map_facet.R"


rule plot_dapc_map_facet:
    input:
        unpack(_dapc_map_facet_inputs),
    output:
        pdf="results/{project}/dapc/plots/{project}.dapc.map-facet.pdf",
        rds="results/{project}/dapc/plots/{project}.dapc.map-facet.rds",
    params:
        lambda wildcards: _map_facet_plot_params(wildcards, "DAPC"),
    log:
        "logs/{project}/plot_dapc_map_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_map_facet.R"


rule plot_snmf_map_facet:
    input:
        unpack(_snmf_map_facet_inputs),
    output:
        pdf="results/{project}/snmf/plots/{project}.snmf.map-facet.pdf",
        rds="results/{project}/snmf/plots/{project}.snmf.map-facet.rds",
    params:
        lambda wildcards: _map_facet_plot_params(wildcards, "sNMF"),
    log:
        "logs/{project}/plot_snmf_map_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_map_facet.R"


rule plot_tess3_map_facet:
    input:
        unpack(_tess3_map_facet_inputs),
    output:
        pdf="results/{project}/tess3/plots/{project}.tess3.map-facet.pdf",
        rds="results/{project}/tess3/plots/{project}.tess3.map-facet.rds",
    params:
        lambda wildcards: _map_facet_plot_params(wildcards, "tess3"),
    log:
        "logs/{project}/plot_tess3_map_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_map_facet.R"


rule plot_alstructure_map_facet:
    input:
        unpack(_alstructure_map_facet_inputs),
    output:
        pdf="results/{project}/alstructure/plots/{project}.alstructure.map-facet.pdf",
        rds="results/{project}/alstructure/plots/{project}.alstructure.map-facet.rds",
    params:
        lambda wildcards: _map_facet_plot_params(wildcards, "ALStructure"),
    log:
        "logs/{project}/plot_alstructure_map_facet.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_map_facet.R"
