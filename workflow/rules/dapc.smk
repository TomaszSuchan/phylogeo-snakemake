# Rule to extract BIC values from find.clusters
# Extracts BIC values from Kstat and writes them to log file
# Uses a shell wrapper with timeout to ensure log is written even if script hangs at prompt
rule dapc_bic_plot:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        log_file = "results/{project}/dapc/{project}.dapc.bic_plot.log.txt"
    params:
        k_values = lambda wildcards: config["projects"][wildcards.project]["parameters"]["k_values"],
        k_values_str = lambda wildcards: ",".join(map(str, config["projects"][wildcards.project]["parameters"]["k_values"])),
        n_pca = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("dapc", {}).get("n_pca", "retained"),
        criterion = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("dapc", {}).get("criterion", "diffNgroup"),
        script = "scripts/dapc_bic_plot.R"
    conda:
        "../envs/dapc.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["dapc"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["dapc"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["dapc"]["runtime"]
    log:
        "logs/{project}/dapc_bic_plot.log"
    shell:
        """
        # Use Python wrapper script with timeout (works on macOS and Linux)
        # The BIC values are extracted and written to log file, no PDF is generated
        # Snakemake runs from project root, so use workflow/scripts/ path
        # Pass arguments to R script: vcf_file, log_file, k_values, n_pca, criterion
        # Use pre-formatted k_values string from params
        python3 workflow/scripts/run_with_timeout.py 90 {output.log_file} workflow/{params.script} {input.vcf} {output.log_file} {params.k_values_str} {params.n_pca} {params.criterion} 2>&1 | tee {log}
        """

# DAPC (Discriminant Analysis of Principal Components) rule - run for each K
rule dapc_analysis:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        results_rds = "results/{project}/dapc/{project}.dapc.K{k}.results.rds",
        cluster_assignments = "results/{project}/dapc/{project}.dapc.K{k}.cluster_assignments.txt",
        membership_probs = "results/{project}/dapc/{project}.dapc.K{k}.membership_probs.txt",
        scatter_plot = "results/{project}/dapc/plots/{project}.dapc.K{k}.scatter.pdf",
        scatter_plot_rds = "results/{project}/dapc/plots/{project}.dapc.K{k}.scatter.rds",
        log_file = "results/{project}/dapc/{project}.dapc.K{k}.log.txt"
    params:
        k = lambda wildcards: int(wildcards.k),
        n_pca = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("dapc", {}).get("n_pca", "retained"),
        n_da = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("dapc", {}).get("n_da", "all"),
        criterion = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("dapc", {}).get("criterion", "diffNgroup")
    conda:
        "../envs/dapc.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["dapc"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["dapc"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["dapc"]["runtime"]
    log:
        "logs/{project}/dapc.K{k}.log"
    benchmark:
        "benchmarks/{project}/dapc.K{k}.txt"
    script:
        "../scripts/dapc_analysis.R"

# Rule to plot BIC values from native find.clusters output
rule dapc_bic_plot_from_log:
    input:
        log_file = rules.dapc_bic_plot.output.log_file
    output:
        plot = "results/{project}/dapc/plots/{project}.dapc.criterion_plot.pdf",
        plot_rds = "results/{project}/dapc/plots/{project}.dapc.criterion_plot.rds"
    conda:
        "../envs/dapc.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    log:
        "logs/{project}/dapc_bic_plot_from_log.log"
    script:
        "../scripts/plot_dapc_bic_from_log.R"

# Rule to plot DAPC results on map with mapmixture (for specific K)
rule mapmixture_dapc:
    input:
        dapc_results = rules.dapc_analysis.output.results_rds,
        popmap = rules.generate_popmap.output.popmap,
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_mapmixture.output
    output:
        plot = "results/{project}/dapc/plots/{project}.dapc.K{k}.map.pdf",
        plot_rds = "results/{project}/dapc/plots/{project}.dapc.K{k}.map.rds"
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
        # Mapmixture-specific parameters
        pie_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_size", 1),
        pie_border = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border", 0.2),
        pie_border_col = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_border_col", "black"),
        pie_opacity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("pie_opacity", 1),
        legend = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("legend", False),
        structure_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"]["mapmixture"].get("structure_colors", ["#66CCEE", "#EE6677", "#228833", "#CCBB44", "#AA3377", "#4477AA", "#BBBBBB", "#EE9988", "#88CCEE", "#CC6677"])
    log:
        "logs/{project}/mapmixture_dapc.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_dapc_map.R"
