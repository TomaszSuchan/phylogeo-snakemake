# Facet panels combining per-K ancestry barplot RDS objects (K >= 2).

def _barplot_facet_plot_params(wildcards, method_label):
    """Layout for multi-K barplot facet panels (single column)."""
    mm = config["projects"][wildcards.project]["parameters"]["mapmixture"]
    facet = config["projects"][wildcards.project]["parameters"]["barplot_facet_plot"]
    return {
        "method_label": method_label,
        "label_width": _fig_cm_to_in(facet["label_width"]),
        "panel_gap": facet["panel_gap"],
        "legend_pad": _fig_cm_to_in(facet.get("legend_pad"), 1.016),
        "flip_axis": mm.get("flip_axis", False),
    }


def _barplot_facet_vertical_plot_params(wildcards, method_label):
    """Layout for multi-K *vertical* barplot facet panels (single row).

    Panels are drawn as vertical barplots (individuals top-to-bottom) tiled
    left-to-right with the lowest K on the left and the highest K on the right.
    Shares label_width/panel_gap/legend_pad with the column layout and adds two
    optional knobs (col_width in cm, label_size in pt) with script-side defaults.
    """
    facet = config["projects"][wildcards.project]["parameters"]["barplot_facet_plot"]
    col_width_cm = facet.get("vertical_col_width")
    # The column-layout label_width sizes a narrow rot-90 "K = n" strip, which is
    # far too small for horizontal site-name labels / the cluster key here, so use
    # generous vertical-specific defaults (still overridable via config).
    return {
        "method_label": method_label,
        "label_width": _fig_cm_to_in(facet.get("vertical_label_width"), 3.8),
        "panel_gap": facet["panel_gap"],
        "legend_pad": _fig_cm_to_in(facet.get("vertical_legend_pad"), 2.5),
        "col_width": _fig_cm_to_in(col_width_cm) if col_width_cm is not None else None,
        "label_size": facet.get("vertical_label_size"),
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


# ---------------------------------------------------------------------------
# Vertical facet panels: per-K barplots tiled left-to-right (lowest K on the
# left, highest K on the right), each drawn as a vertical barplot. Reuse the
# same per-K barplot RDS inputs as the column layout above.
# ---------------------------------------------------------------------------


rule plot_structure_barplot_facet_vertical:
    input:
        unpack(_structure_barplot_facet_inputs),
    output:
        pdf="results/{project}/structure/plots/{project}.structure.barplot-facet-vertical.pdf",
        rds="results/{project}/structure/plots/{project}.structure.barplot-facet-vertical.rds",
    params:
        lambda wildcards: _barplot_facet_vertical_plot_params(wildcards, "STRUCTURE"),
    log:
        "logs/{project}/plot_structure_barplot_facet_vertical.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet_vertical.R"


rule plot_faststructure_barplot_facet_vertical:
    input:
        unpack(_faststructure_barplot_facet_inputs),
    output:
        pdf="results/{project}/faststructure/plots/{project}.faststructure.barplot-facet-vertical.pdf",
        rds="results/{project}/faststructure/plots/{project}.faststructure.barplot-facet-vertical.rds",
    params:
        lambda wildcards: _barplot_facet_vertical_plot_params(wildcards, "fastStructure"),
    log:
        "logs/{project}/plot_faststructure_barplot_facet_vertical.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet_vertical.R"


rule plot_admixture_barplot_facet_vertical:
    input:
        unpack(_admixture_barplot_facet_inputs),
    output:
        pdf="results/{project}/admixture/plots/{project}.admixture.barplot-facet-vertical.pdf",
        rds="results/{project}/admixture/plots/{project}.admixture.barplot-facet-vertical.rds",
    params:
        lambda wildcards: _barplot_facet_vertical_plot_params(wildcards, "ADMIXTURE"),
    log:
        "logs/{project}/plot_admixture_barplot_facet_vertical.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet_vertical.R"


rule plot_snmf_barplot_facet_vertical:
    input:
        unpack(_snmf_barplot_facet_inputs),
    output:
        pdf="results/{project}/snmf/plots/{project}.snmf.barplot-facet-vertical.pdf",
        rds="results/{project}/snmf/plots/{project}.snmf.barplot-facet-vertical.rds",
    params:
        lambda wildcards: _barplot_facet_vertical_plot_params(wildcards, "sNMF"),
    log:
        "logs/{project}/plot_snmf_barplot_facet_vertical.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet_vertical.R"


rule plot_tess3_barplot_facet_vertical:
    input:
        unpack(_tess3_barplot_facet_inputs),
    output:
        pdf="results/{project}/tess3/plots/{project}.tess3.barplot-facet-vertical.pdf",
        rds="results/{project}/tess3/plots/{project}.tess3.barplot-facet-vertical.rds",
    params:
        lambda wildcards: _barplot_facet_vertical_plot_params(wildcards, "tess3"),
    log:
        "logs/{project}/plot_tess3_barplot_facet_vertical.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet_vertical.R"


rule plot_alstructure_barplot_facet_vertical:
    input:
        unpack(_alstructure_barplot_facet_inputs),
    output:
        pdf="results/{project}/alstructure/plots/{project}.alstructure.barplot-facet-vertical.pdf",
        rds="results/{project}/alstructure/plots/{project}.alstructure.barplot-facet-vertical.rds",
    params:
        lambda wildcards: _barplot_facet_vertical_plot_params(wildcards, "ALStructure"),
    log:
        "logs/{project}/plot_alstructure_barplot_facet_vertical.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet_vertical.R"


rule plot_construct_barplot_facet_vertical:
    input:
        unpack(_construct_barplot_facet_inputs),
    output:
        pdf="results/{project}/construct/plots/{project}.construct.barplot-facet-vertical.pdf",
        rds="results/{project}/construct/plots/{project}.construct.barplot-facet-vertical.rds",
    params:
        lambda wildcards: _barplot_facet_vertical_plot_params(wildcards, "conStruct"),
    log:
        "logs/{project}/plot_construct_barplot_facet_vertical.log",
    conda:
        "../envs/r-plot.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_structure_barplot_facet_vertical.R"
