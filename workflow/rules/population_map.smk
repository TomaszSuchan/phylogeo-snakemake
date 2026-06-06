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


# DEM basemap via elevatr (cached in EPSG:4326 for mapmixture cropping).
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
        plot = "results/{project}/stats_samples/plots/{project}.population_map.pdf",
        plot_rds = "results/{project}/stats_samples/plots/{project}.population_map.rds"
    params:
        lambda wildcards: _population_map_params(wildcards)
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

