rule prepare_invariant_vcf:
    input:
        loci = lambda wildcards: config["projects"][wildcards.project]["ipyrad_prefix"] + ".loci",
        samples_file = rules.update_samples_file.output.filtered_samples
    output:
        invariant_vcf = "results/{project}/filtered_data/{project}.invariant_sites.vcf"
    log:
        "logs/{project}/prepare_invariant_vcf.log"
    benchmark:
        "benchmarks/{project}/prepare_invariant_vcf.txt"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/extract_invariant_vcf.py {input.loci} -o {output.invariant_vcf} --samples-file {input.samples_file} &> {log}
        """

rule prepare_invariant_vcf_gz:
    input:
        vcf = rules.prepare_invariant_vcf.output.invariant_vcf
    output:
        vcf = temporary("results/{project}/filtered_data/{project}.invariant_sites.vcf.gz")
    log:
        "logs/{project}/prepare_invariant_vcf_gz.log"
    benchmark:
        "benchmarks/{project}/prepare_invariant_vcf_gz.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bgzip {input.vcf} &> {log}
        """

rule prepare_invariant_vcf_gz_index:
    input:
        vcf = rules.prepare_invariant_vcf_gz.output.vcf
    output:
        index = temporary("results/{project}/filtered_data/{project}.invariant_sites.vcf.gz.csi")
    log:
        "logs/{project}/prepare_invariant_vcf_gz_index.log"
    benchmark:
        "benchmarks/{project}/prepare_invariant_vcf_gz_index.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools index {input.vcf}  &> {log}
        """

rule merge_invariant_sites:
    input:
        vcf_var = rules.subset_vcf_after_relatedness.output.vcf,
        vcf_var_index = rules.index_vcf_after_relatedness.output.index,
        vcf_inv = rules.prepare_invariant_vcf_gz.output.vcf,
        vcf_inv_index = rules.prepare_invariant_vcf_gz_index.output.index
    output:
        merged_vcf = "results/{project}/filtered_data/{project}.merged_invariant_sites.vcf.gz",
        index = "results/{project}/filtered_data/{project}.merged_invariant_sites.vcf.gz.csi"
    log:
        "logs/{project}/merge_invariant_sites.log"
    benchmark:
        "benchmarks/{project}/merge_invariant_sites.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools concat {input.vcf_inv} {input.vcf_var} -Oz -a -o {output.merged_vcf}  >> {log} 2>&1
        bcftools index {output.merged_vcf}
        """

rule pixy:
    input:
        vcf = rules.merge_invariant_sites.output.merged_vcf,
        popmap = rules.generate_popmap.output
    output:
        pi = "results/{project}/pixy/{project}.pixy_pi.txt",
        fst = "results/{project}/pixy/{project}.pixy_fst.txt",
        dxy = "results/{project}/pixy/{project}.pixy_dxy.txt"
    log:
        "logs/{project}/pixy.log"
    benchmark:
        "benchmarks/{project}/pixy.txt"
    params:
        window_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("window_size", 10000),
        output_folder = "results/{project}/pixy/",
        output_prefix = "{project}.pixy"
    conda:
        "../envs/pixy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["runtime"]
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats pi fst dxy --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder} --output_prefix {params.output_prefix} &> {log}
        """

rule pixy_summary:
    input:
        pi = rules.pixy.output.pi,
        fst = rules.pixy.output.fst,
        dxy = rules.pixy.output.dxy
    output:
        pi = "results/{project}/pixy/{project}.pixy_pi-summary.txt",
        fst = "results/{project}/pixy/{project}.pixy_fst-summary.txt",
        dxy = "results/{project}/pixy/{project}.pixy_dxy-summary.txt"
    log:
        "logs/{project}/pixy_summary.log"
    benchmark:
        "benchmarks/{project}/pixy_summary.txt"
    params:
        bootstrap_replicates = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("bootstrap_replicates", 1000)
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/pixy_summary.py"

# Rule to plot FST heatmap with dendrogram
rule plot_pixy_fst_heatmap:
    input:
        fst_summary = rules.pixy_summary.output.fst
    output:
        pdf = "results/{project}/pixy/plots/{project}.pixy_fst_heatmap.pdf",
        rds = "results/{project}/pixy/plots/{project}.pixy_fst_heatmap.rds"
    log:
        "logs/{project}/plot_pixy_fst_heatmap.log"
    benchmark:
        "benchmarks/{project}/plot_pixy_fst_heatmap.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_pixy_fst_heatmap.R"

# Rule to plot DXY heatmap with dendrogram
rule plot_pixy_dxy_heatmap:
    input:
        dxy_summary = rules.pixy_summary.output.dxy
    output:
        pdf = "results/{project}/pixy/plots/{project}.pixy_dxy_heatmap.pdf",
        rds = "results/{project}/pixy/plots/{project}.pixy_dxy_heatmap.rds"
    log:
        "logs/{project}/plot_pixy_dxy_heatmap.log"
    benchmark:
        "benchmarks/{project}/plot_pixy_dxy_heatmap.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_pixy_dxy_heatmap.R"

# Rule to plot Pi barplot with confidence intervals
# When color_by != "none", generates grouped version (grouped by stratification, then sorted by pi within group)
rule plot_pixy_pi_barplot:
    input:
        pi_summary = rules.pixy_summary.output.pi,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/pixy/plots/{project}.pixy_pi_barplot-grouped-{color_by}.pdf",
        rds = "results/{project}/pixy/plots/{project}.pixy_pi_barplot-grouped-{color_by}.rds"
    log:
        "logs/{project}/plot_pixy_pi_barplot_{color_by}.log"
    params:
        color_by = lambda wildcards: wildcards.color_by if wildcards.color_by != "none" else None,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        plot_type = "grouped"  # Grouped by stratification
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_pixy_pi_barplot.R"

# Rule to plot Pi barplot sorted by pi (colored by stratification but not grouped)
# Only generated when color_by != "none"
rule plot_pixy_pi_barplot_sorted:
    input:
        pi_summary = rules.pixy_summary.output.pi,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/pixy/plots/{project}.pixy_pi_barplot-sorted-{color_by}.pdf",
        rds = "results/{project}/pixy/plots/{project}.pixy_pi_barplot-sorted-{color_by}.rds"
    log:
        "logs/{project}/plot_pixy_pi_barplot_{color_by}_sorted.log"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        plot_type = "sorted"  # Sorted by pi value
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_pixy_pi_barplot.R"

# Rule to plot Pi barplot without coloring (plain barplot)
rule plot_pixy_pi_barplot_plain:
    input:
        pi_summary = rules.pixy_summary.output.pi,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/pixy/plots/{project}.pixy_pi_barplot-plain.pdf",
        rds = "results/{project}/pixy/plots/{project}.pixy_pi_barplot-plain.rds"
    log:
        "logs/{project}/plot_pixy_pi_barplot_plain.log"
    params:
        color_by = None,  # No coloring
        pca_colors = None,
        plot_type = "plain"  # Plain barplot, sorted by pi value
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_pixy_pi_barplot.R"

# Rule to plot Pi on map
rule plot_pixy_pi_map:
    input:
        popmap = rules.generate_popmap.output.popmap,
        indpopdata = rules.generate_popdata.output.indpopdata,
        summary = rules.pixy_summary.output.pi,
        install = rules.install_mapmixture.output  # Reuse mapmixture installation
    output:
        pdf = "results/{project}/pixy/plots/{project}.pixy_pi_map.pdf",
        rds = "results/{project}/pixy/plots/{project}.pixy_pi_map.rds"
    params:
        stat_type = "pi",
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
        # Pixy-specific parameters
        point_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("point_size", 3)
    log:
        "logs/{project}/plot_pixy_pi_map.log"
    benchmark:
        "benchmarks/{project}/plot_pixy_pi_map.txt"
    conda:
        "../envs/mapmixture.yaml"  # Reuse mapmixture environment
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_pixy_maps.R"
