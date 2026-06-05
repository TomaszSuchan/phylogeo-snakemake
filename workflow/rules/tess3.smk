# tess3r spatial ancestry inference and plotting


def _tess3_params(wildcards):
    return config["projects"][wildcards.project]["parameters"].get("tess3", {})


rule install_tess3:
    """
    Install tess3r once into the Snakemake conda environment.
    tess3r is installed from its upstream GitHub repository after conda
    creates the base R/vcfR plotting environment.
    """
    output:
        touch(".snakemake/tess3_installed")
    conda:
        "../envs/tess3.yaml"
    shell:
        """
        Rscript --vanilla -e 'lib <- .libPaths()[1]; unlink(list.files(lib, pattern="^00LOCK-", full.names=TRUE), recursive=TRUE, force=TRUE); if (!requireNamespace("tess3r", quietly=TRUE, lib.loc=lib)) remotes::install_github("bcm-uga/TESS3_encho_sen", upgrade="never", dependencies=TRUE, lib=lib); if (!requireNamespace("tess3r", quietly=TRUE, lib.loc=lib)) stop("Failed to install required R package: tess3r"); options(repos="https://cloud.r-project.org/"); pkgs <- c("mapmixture", "elevatr"); for (p in pkgs) if (!requireNamespace(p, quietly=TRUE, lib.loc=lib)) install.packages(p, lib=lib)'
        """


rule tess3_genotypes:
    """
    Convert the analysis VCF to a tess3r-compatible genotype dosage matrix.
    """
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards),
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_tess3.output
    output:
        geno_rds = temporary("results/{project}/tess3/{project}.genotypes.rds")
    log:
        "logs/{project}/tess3_genotypes.log"
    benchmark:
        "benchmarks/{project}/tess3_genotypes.txt"
    params:
        ploidy = lambda wildcards: _tess3_params(wildcards).get("ploidy", 2)
    conda:
        "../envs/tess3.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["runtime"]
    script:
        "../scripts/tess3_genotypes.R"


rule tess3_choose_k:
    """
    Run tess3r across all configured K values with masked cross-validation.
    """
    input:
        geno_rds = rules.tess3_genotypes.output.geno_rds,
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_tess3.output
    output:
        results_rds = "results/{project}/tess3/{project}.tess3.chooseK.results.rds",
        choose_k_results = "results/{project}/tess3/{project}.tess3.chooseK_results.txt",
        cv_summary = "results/{project}/tess3/{project}.tess3.cv_summary.txt"
    log:
        "logs/{project}/tess3_choose_k.log"
    benchmark:
        "benchmarks/{project}/tess3_choose_k.txt"
    params:
        k_values = lambda wildcards: config["projects"][wildcards.project]["parameters"]["k_values"],
        method = lambda wildcards: _tess3_params(wildcards).get("method", "projected.ls"),
        replicates = lambda wildcards: _tess3_params(wildcards).get("replicates", 1),
        max_iteration = lambda wildcards: _tess3_params(wildcards).get("max_iteration", 200),
        tolerance = lambda wildcards: _tess3_params(wildcards).get("tolerance", 1e-05),
        ploidy = lambda wildcards: _tess3_params(wildcards).get("ploidy", 2),
        mask = lambda wildcards: _tess3_params(wildcards).get("mask", 0),
        crossvalid = lambda wildcards: _tess3_params(wildcards).get("crossvalid", False),
        crossentropy = lambda wildcards: _tess3_params(wildcards).get("crossentropy", False)
    conda:
        "../envs/tess3.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["runtime"]
    script:
        "../scripts/tess3_choose_k.R"


rule plot_tess3_cv:
    """
    Plot tess3r cross-validation scores across K (evanno-style ggplot aesthetics).
    """
    input:
        cv_summary = rules.tess3_choose_k.output.cv_summary
    output:
        pdf = "results/{project}/tess3/plots/{project}.tess3.cv_plot.pdf",
        rds = "results/{project}/tess3/plots/{project}.tess3.cv_plot.rds"
    log:
        "logs/{project}/plot_tess3_cv.log"
    params:
        crossvalid = lambda wildcards: _tess3_params(wildcards).get("crossvalid", False),
        crossentropy = lambda wildcards: _tess3_params(wildcards).get("crossentropy", False),
        width = lambda wildcards: _choose_k_plot_param(wildcards, "width", 10),
        height = lambda wildcards: _choose_k_plot_param(wildcards, "height", 5),
        dpi = lambda wildcards: _choose_k_plot_param(wildcards, "dpi", 300)
    conda:
        "../envs/tess3.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_choose_k_score.R"


rule tess3_analysis:
    """
    Run tess3r spatially regularized ancestry estimation for a given K.
    """
    input:
        geno_rds = rules.tess3_genotypes.output.geno_rds,
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
        replicates = lambda wildcards: _tess3_params(wildcards).get("replicates", 1),
        max_iteration = lambda wildcards: _tess3_params(wildcards).get("max_iteration", 200),
        tolerance = lambda wildcards: _tess3_params(wildcards).get("tolerance", 1e-05),
        ploidy = lambda wildcards: _tess3_params(wildcards).get("ploidy", 2),
        map_method = lambda wildcards: _tess3_params(wildcards).get("map_method", "map.max"),
        map_resolution = lambda wildcards: _tess3_params(wildcards).get("map_resolution", [300, 300]),
        interpolation_knots = lambda wildcards: _tess3_params(wildcards).get("interpolation_knots", 10)
    conda:
        "../envs/tess3.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["tess3"]["runtime"]
    script:
        "../scripts/tess3_analysis.R"


rule plot_tess3_ggmap:
    """
    Interpolated tess3r ancestry map (ggtess3Q) on the shared mapmixture basemap.
    """
    input:
        unpack(_tess3_ggmap_inputs),
    output:
        pdf = "results/{project}/tess3/plots/{project}.tess3.K{k}.ggmap.pdf",
        rds = "results/{project}/tess3/plots/{project}.tess3.K{k}.ggmap.rds"
    log:
        "logs/{project}/plot_tess3_ggmap.K{k}.log"
    params:
        k = lambda wildcards: int(wildcards.k),
        map_resolution = lambda wildcards: _tess3_params(wildcards).get("map_resolution", [300, 300]),
        interpolation_knots = lambda wildcards: _tess3_params(wildcards).get("interpolation_knots", 10),
        structure_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("structure_colors", ["#66CCEE", "#EE6677", "#228833", "#CCBB44", "#AA3377", "#4477AA", "#BBBBBB", "#EE9988", "#88CCEE", "#CC6677"]),
        point_size = lambda wildcards: _tess3_params(wildcards).get("ggmap_point_size", 0.8),
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
        use_elevation_bg = lambda wildcards: _use_elevation_bg(wildcards)
    conda:
        "../envs/tess3.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_tess3_ggmap.R"


rule plot_tess3_interpolation:
    """
    Interpolated ancestry map using tess3r native plotting (FieldsKrigModel).
    """
    input:
        results_rds = rules.tess3_analysis.output.results_rds
    output:
        pdf = "results/{project}/tess3/plots/{project}.tess3.K{k}.interpolation.pdf"
    log:
        "logs/{project}/plot_tess3_interpolation.K{k}.log"
    params:
        k = lambda wildcards: int(wildcards.k),
        map_method = lambda wildcards: _tess3_params(wildcards).get("map_method", "map.max"),
        map_resolution = lambda wildcards: _tess3_params(wildcards).get("map_resolution", [300, 300]),
        interpolation_knots = lambda wildcards: _tess3_params(wildcards).get("interpolation_knots", 10),
        structure_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("structure_colors", ["#66CCEE", "#EE6677", "#228833", "#CCBB44", "#AA3377", "#4477AA", "#BBBBBB", "#EE9988", "#88CCEE", "#CC6677"])
    conda:
        "../envs/tess3.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_tess3_interpolation.R"


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
    benchmark:
        "benchmarks/{project}/mapmixture_tess3.K{k}.txt"
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
    benchmark:
        "benchmarks/{project}/tess3_barplot.K{k}.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_barplots.R"
