# tess3r spatial ancestry inference and plotting


def _tess3_params(wildcards):
    return config["projects"][wildcards.project]["parameters"].get("tess3", {})


rule install_tess3:
    """
    Install tess3r once into the Snakemake conda environment.
    tess3r is CRAN-only for this workflow, so it is installed after conda
    creates the base R/vcfR plotting environment.
    """
    output:
        touch(".snakemake/tess3_installed")
    conda:
        "../envs/tess3.yaml"
    shell:
        """
        Rscript --vanilla -e 'lib <- .libPaths()[1]; unlink(list.files(lib, pattern="^00LOCK-", full.names=TRUE), recursive=TRUE, force=TRUE); needed <- "tess3r"; missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing) > 0) install.packages(missing, repos="https://cloud.r-project.org"); missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing) > 0) stop("Failed to install required R packages: ", paste(missing, collapse=", "))'
        """


rule tess3_analysis:
    """
    Run tess3r spatially regularized ancestry estimation for a given K.
    """
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards),
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_tess3.output
    output:
        qmatrix = "results/{project}/tess3/{project}.tess3.K{k}.Qmatrix.txt",
        results_rds = "results/{project}/tess3/{project}.tess3.K{k}.results.rds",
        cross_entropy = "results/{project}/tess3/{project}.tess3.K{k}.cross_entropy.txt",
        max_cluster_png = "results/{project}/tess3/plots/{project}.tess3.K{k}.max_cluster.png"
    log:
        "logs/{project}/tess3.K{k}.log"
    benchmark:
        "benchmarks/{project}/tess3.K{k}.txt"
    params:
        k = lambda wildcards: int(wildcards.k),
        method = lambda wildcards: _tess3_params(wildcards).get("method", "projected.ls"),
        replicates = lambda wildcards: _tess3_params(wildcards).get("replicates", 5),
        max_iteration = lambda wildcards: _tess3_params(wildcards).get("max_iteration", 1000),
        tolerance = lambda wildcards: _tess3_params(wildcards).get("tolerance", 1e-05),
        n_colors = lambda wildcards: _tess3_params(wildcards).get("n_colors", 16),
        ploidy = lambda wildcards: _tess3_params(wildcards).get("ploidy", 2)
    conda:
        "../envs/tess3.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["runtime"]
    script:
        "../scripts/tess3_analysis.R"


rule mapmixture_tess3:
    """
    Plot tess3r ancestry coefficients on the shared mapmixture basemap.
    """
    input:
        unpack(_tess3_map_inputs),
    output:
        plot = "results/{project}/tess3/plots/{project}.tess3.K{k}.map.pdf",
        plot_rds = "results/{project}/tess3/plots/{project}.tess3.K{k}.map.rds"
    params:
        output_prefix = "results/{project}/tess3/{project}.tess3.K{k}",
        width = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("width", 10),
        height = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("height", 8),
        dpi = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("dpi", 300),
        boundary = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("map_boundary", "NULL"),
        crs = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("crs", 4326),
        basemap = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("basemap", "NULL"),
        land_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("land_colour", "#d9d9d9"),
        sea_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("sea_colour", "#deebf7"),
        expand = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("expand", False),
        arrow = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("arrow", True),
        arrow_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("arrow_size", 1),
        arrow_position = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("arrow_position", "tl"),
        scalebar = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("scalebar", True),
        scalebar_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("scalebar_size", 1),
        scalebar_position = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("scalebar_position", "tl"),
        plot_title = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("plot_title", ""),
        axis_title_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("axis_title_size", 10),
        axis_text_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("axis_text_size", 8),
        basemap_border = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("basemap_border", True),
        basemap_border_col = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("basemap_border_col", "black"),
        basemap_border_lwd = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("basemap_border_lwd", 0.1),
        use_elevation_bg = lambda wildcards: _use_elevation_bg(wildcards),
        pie_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_size", 1),
        pie_border = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border", 0.2),
        pie_border_col = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border_col", "black"),
        pie_opacity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_opacity", 1),
        legend = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("legend", False),
        structure_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("structure_colors", ["#66CCEE", "#EE6677", "#228833", "#CCBB44", "#AA3377", "#4477AA", "#BBBBBB", "#EE9988", "#88CCEE", "#CC6677"])
    log:
        "logs/{project}/mapmixture_tess3.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"


rule barplot_tess3:
    """
    Draw tess3r ancestry barplots with the shared mapmixture barplot script.
    """
    input:
        qmatrix = rules.tess3_analysis.output.qmatrix,
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = ".snakemake/mapmixture_installed"
    output:
        barplot = "results/{project}/tess3/plots/{project}.tess3.K{k}.barplot.pdf",
        barplot_rds = "results/{project}/tess3/plots/{project}.tess3.K{k}.barplot.rds"
    params:
        output_prefix = "results/{project}/tess3/{project}.tess3.K{k}",
        width = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("width", 10),
        height = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("height", 6),
        dpi = lambda wildcards: config["projects"][wildcards.project]["parameters"]["map_background"].get("dpi", 300),
        site_dividers = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("site_dividers", True),
        divider_width = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("divider_width", 0.4),
        site_order = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("site_order", "NULL"),
        population_sort_by = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("population_sort_by", "NULL"),
        flip_axis = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("flip_axis", False),
        site_labels_angle = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("site_labels_angle", 90),
        population_labels = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("population_labels", ["Site"]),
        structure_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("structure_colors", ["#66CCEE", "#EE6677", "#228833", "#CCBB44", "#AA3377", "#4477AA", "#BBBBBB", "#EE9988", "#88CCEE", "#CC6677"])
    log:
        "logs/{project}/tess3_barplot.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_barplots.R"
