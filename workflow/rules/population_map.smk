# Basemap caches for population maps (and shared map outputs under results/{project}/maps/).
# Requires install_mapmixture from mapmixture.smk — include mapmixture.smk before this file in the Snakefile.

# Download + extract global Natural Earth 10m HR zip once (params: ne_kind only — not CRS/bbox).
# ne1: NE1_HR_LC_SR_W_DR.zip; gray_earth: GRAY_HR_SR_OB_DR.zip (naciscdn.org naturalearth/10m/raster).
rule fetch_naturalearth_hr_global:
    input:
        install=rules.install_mapmixture.output,
    output:
        zip="results/{project}/maps/naturalearth_hr_global/hr_bundle.zip",
        tif="results/{project}/maps/naturalearth_hr_global/hr_global.tif",
    params:
        ne_kind=lambda w: _natural_earth_kind(w),
    log:
        "logs/{project}/fetch_naturalearth_hr_global.log",
    conda:
        "../envs/mapmixture.yaml"
    threads: 1
    resources:
        mem_mb=lambda w: config["projects"][w.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda w: max(
            180,
            int(config["projects"][w.project]["parameters"]["resources"]["default"]["runtime"]),
        ),
    script:
        "../scripts/fetch_naturalearth_hr_global.R"


# Crop + warp global HR tif to study extent and map_background.crs (reruns when CRS/bbox/indpopdata change).
rule cache_naturalearth_basemap:
    input:
        indpopdata=rules.generate_popdata.output.indpopdata,
        install=rules.install_mapmixture.output,
        hr_global=rules.fetch_naturalearth_hr_global.output.tif,
    output:
        "results/{project}/maps/naturalearth_basemap.tif",
    params:
        boundary=lambda w: config["projects"][w.project]["parameters"].get("map_boundary", "NULL"),
        crs=lambda w: config["projects"][w.project]["parameters"]["map_background"].get("crs", 4326),
    log:
        "logs/{project}/cache_naturalearth_basemap.log",
    conda:
        "../envs/mapmixture.yaml"
    threads: 1
    resources:
        mem_mb=lambda w: config["projects"][w.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda w: config["projects"][w.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/cache_naturalearth_basemap.R"


# DEM basemap via elevatr (written in map_background.crs).
rule cache_map_elevation:
    input:
        indpopdata=rules.generate_popdata.output.indpopdata,
        install=rules.install_mapmixture.output,
    output:
        "results/{project}/maps/elevation_basemap.tif",
    params:
        boundary=lambda w: config["projects"][w.project]["parameters"].get("map_boundary", "NULL"),
        crs=lambda w: config["projects"][w.project]["parameters"]["map_background"].get("crs", 4326),
        width=lambda w: config["projects"][w.project]["parameters"]["map_background"].get("width", 10),
        height=lambda w: config["projects"][w.project]["parameters"]["map_background"].get("height", 8),
        dpi=lambda w: config["projects"][w.project]["parameters"]["map_background"].get("dpi", 300),
        elevatr_z=lambda w: config["projects"][w.project]["parameters"]["map_background"].get("elevatr_z"),
    log:
        "logs/{project}/cache_map_elevation.log",
    conda:
        "../envs/mapmixture.yaml"
    threads: 1
    resources:
        mem_mb=lambda w: config["projects"][w.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda w: config["projects"][w.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/cache_map_elevation.R"


# Rule for plotting population locations on map with labels and leader lines

rule plot_population_map:
    input:
        unpack(_population_map_inputs),
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
        use_elevation_bg=lambda wildcards: _use_elevation_bg(wildcards),
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

