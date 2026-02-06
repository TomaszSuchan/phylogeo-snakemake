# Rule for plotting population locations on map with labels and leader lines

rule plot_population_map:
    input:
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_mapmixture.output  # Reuse mapmixture installation
    output:
        plot = "results/{project}/stats_samples/{project}.population_map.pdf"
    params:
        # Map background parameters (shared)
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
        # Label-specific parameters
        point_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("point_size", 3),
        point_color = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("point_color", "black"),
        point_shape = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("point_shape", 19),
        label_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("label_size", 3),
        label_color = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("label_color", "black"),
        label_fontface = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("label_fontface", "plain"),
        show_points = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("show_points", True),
        show_labels = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("show_labels", True),
        # ggrepel parameters
        force = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("force", 10),
        force_pull = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("force_pull", 1),
        max_overlaps = lambda wildcards: str(config["projects"][wildcards.project]["parameters"]["population_map"].get("max.overlaps", float("inf"))),
        min_segment_length = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("min.segment.length", 0.5),
        segment_color = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("segment.color", "grey50"),
        segment_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["population_map"].get("segment.size", 0.5)
    log:
        "logs/{project}/plot_population_map.log"
    benchmark:
        "benchmarks/{project}/plot_population_map.txt"
    conda:
        "../envs/mapmixture.yaml"  # Reuse mapmixture environment
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_population_map.R"

