# Rules for plotting STRUCTURE, fastStructure, and ADMIXTURE results on maps using mapmixture

# Rule to install mapmixture once (all plotting rules depend on this)
rule install_mapmixture:
    output:
        touch(".snakemake/mapmixture_installed")
    conda:
        "../envs/mapmixture.yaml"
    shell:
        """
        Rscript -e 'if (!require("mapmixture", quietly = TRUE)) install.packages("mapmixture", repos = "https://cloud.r-project.org/")'
        """

# Rule to plot STRUCTURE results on map with mapmixture
rule mapmixture_structure:
    input:
        qmatrix = "results/{project}/structure/{project}.structure.K{k}.Qmatrix.txt",
        popmap = rules.generate_popmap.output.popmap,
        install = rules.install_mapmixture.output
    output:
        plot = "results/{project}/structure/{project}.structure.K{k}.map.pdf"
    params:
        popdata=lambda wildcards: config["projects"][wildcards.project]["parameters"]["popdata"],
        output_prefix = "results/{project}/structure/{project}.structure.K{k}",
        width = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("width", 10),
        height = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("height", 8),
        dpi = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("dpi", 300),
        boundary = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("boundary", "NULL"),
        crs = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("crs", 4326),
        basemap = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap", "NULL"),
        pie_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_size", 1),
        pie_border = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border", 0.2),
        pie_border_col = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border_col", "black"),
        pie_opacity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_opacity", 1),
        land_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("land_colour", "#d9d9d9"),
        sea_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("sea_colour", "#deebf7"),
        expand = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("expand", False),
        arrow = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("arrow", True),
        arrow_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("arrow_size", 1),
        arrow_position = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("arrow_position", "tl"),
        scalebar = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("scalebar", True),
        scalebar_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("scalebar_size", 1),
        scalebar_position = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("scalebar_position", "tl"),
        plot_title = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("plot_title", ""),
        axis_title_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("axis_title_size", 10),
        axis_text_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("axis_text_size", 8),
        basemap_border = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap_border", True),
        basemap_border_col = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap_border_col", "black"),
        basemap_border_lwd = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap_border_lwd", 0.1),
        legend = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("legend", False)
    log:
        "logs/{project}/mapmixture_structure.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"

# Rule to plot fastStructure results on map with mapmixture
rule mapmixture_faststructure:
    input:
        qmatrix = "results/{project}/faststructure/{project}.faststructure.{k}.meanQ",
        popmap = rules.generate_popmap.output.popmap,
        install = rules.install_mapmixture.output
    output:
        plot = "results/{project}/faststructure/{project}.faststructure.K{k}.map.pdf"
    params:
        popdata=lambda wildcards: config["projects"][wildcards.project]["parameters"]["popdata"],
        output_prefix = "results/{project}/faststructure/{project}.faststructure.K{k}",
        width = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("width", 10),
        height = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("height", 8),
        dpi = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("dpi", 300),
        boundary = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("boundary", "NULL"),
        crs = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("crs", 4326),
        basemap = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap", "NULL"),
        pie_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_size", 1),
        pie_border = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border", 0.2),
        pie_border_col = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border_col", "black"),
        pie_opacity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_opacity", 1),
        land_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("land_colour", "#d9d9d9"),
        sea_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("sea_colour", "#deebf7"),
        expand = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("expand", False),
        arrow = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("arrow", True),
        arrow_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("arrow_size", 1),
        arrow_position = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("arrow_position", "tl"),
        scalebar = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("scalebar", True),
        scalebar_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("scalebar_size", 1),
        scalebar_position = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("scalebar_position", "tl"),
        plot_title = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("plot_title", ""),
        axis_title_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("axis_title_size", 10),
        axis_text_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("axis_text_size", 8),
        basemap_border = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap_border", True),
        basemap_border_col = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap_border_col", "black"),
        basemap_border_lwd = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap_border_lwd", 0.1),
        legend = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("legend", False)
    log:
        "logs/{project}/mapmixture_faststructure.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"

# Rule to plot ADMIXTURE results on map with mapmixture
rule mapmixture_admixture:
    input:
        qmatrix = "results/{project}/admixture/{project}.biallelic_snps_thinned.{k}.Q",
        popmap = rules.generate_popmap.output.popmap,
        install = rules.install_mapmixture.output
    output:
        plot = "results/{project}/admixture/{project}.admixture.K{k}.map.pdf"
    params:
        popdata=lambda wildcards: config["projects"][wildcards.project]["parameters"]["popdata"],
        output_prefix = "results/{project}/admixture/{project}.admixture.K{k}",
        width = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("width", 10),
        height = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("height", 8),
        dpi = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("dpi", 300),
        boundary = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("boundary", "NULL"),
        crs = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("crs", 4326),
        basemap = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap", "NULL"),
        pie_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_size", 1),
        pie_border = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border", 0.2),
        pie_border_col = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border_col", "black"),
        pie_opacity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_opacity", 1),
        land_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("land_colour", "#d9d9d9"),
        sea_colour = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("sea_colour", "#deebf7"),
        expand = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("expand", False),
        arrow = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("arrow", True),
        arrow_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("arrow_size", 1),
        arrow_position = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("arrow_position", "tl"),
        scalebar = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("scalebar", True),
        scalebar_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("scalebar_size", 1),
        scalebar_position = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("scalebar_position", "tl"),
        plot_title = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("plot_title", ""),
        axis_title_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("axis_title_size", 10),
        axis_text_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("axis_text_size", 8),
        basemap_border = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap_border", True),
        basemap_border_col = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap_border_col", "black"),
        basemap_border_lwd = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("basemap_border_lwd", 0.1),
        legend = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("legend", False)
    log:
        "logs/{project}/mapmixture_admixture.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"
