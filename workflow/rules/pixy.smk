rule prepare_invariant_vcf:
    input:
        loci = lambda wildcards: config["projects"][wildcards.project]["ipyrad_prefix"] + ".loci",
        template_vcf = rules.sort_vcf.output.vcf,
        # Use canonical keep list (already includes relatedness filtering when enabled)
        samples_file = rules.filter_related_individuals.output.samples_to_keep
    output:
        invariant_vcf = temporary("results/{project}/filtered_data/{project}.invariant_sites.vcf")
    log:
        "logs/{project}/prepare_invariant_vcf.log"
    benchmark:
        "benchmarks/{project}/prepare_invariant_vcf.txt"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["prepare_invariant_vcf"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["prepare_invariant_vcf"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["prepare_invariant_vcf"]["runtime"]
    shell:
        """
        python workflow/scripts/extract_invariant_vcf.py {input.loci} -o {output.invariant_vcf} --samples-file {input.samples_file} --template-vcf {input.template_vcf} &> {log}
        """

rule prepare_invariant_vcf_gz:
    input:
        vcf = rules.prepare_invariant_vcf.output.invariant_vcf
    output:
        vcf = "results/{project}/filtered_data/{project}.invariant_sites.vcf.gz"
    log:
        "logs/{project}/prepare_invariant_vcf_gz.log"
    benchmark:
        "benchmarks/{project}/prepare_invariant_vcf_gz.txt"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        vcf-sort {input.vcf} 2> {log} | bgzip -c > {output.vcf}
        """

rule prepare_invariant_vcf_gz_index:
    input:
        vcf = rules.prepare_invariant_vcf_gz.output.vcf
    output:
        index = "results/{project}/filtered_data/{project}.invariant_sites.vcf.gz.csi"
    log:
        "logs/{project}/prepare_invariant_vcf_gz_index.log"
    benchmark:
        "benchmarks/{project}/prepare_invariant_vcf_gz_index.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        bcftools index {input.vcf}  &> {log}
        """

rule generate_pixy_popmap:
    """
    Create a pixy popmap from indpopdata, using a configurable grouping column.
    """
    input:
        indpopdata = rules.generate_popdata.output.indpopdata
    output:
        popmap = "results/{project}/pixy/{project}.{grouping}.pixy_popmap.txt"
    params:
        group_by = lambda wildcards: wildcards.grouping
    log:
        "logs/{project}/generate_pixy_popmap.{grouping}.log"
    benchmark:
        "benchmarks/{project}/generate_pixy_popmap.{grouping}.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/generate_pixy_popmap_from_indpopdata.py"

rule pixy_pi:
    input:
        vcf = rules.prepare_invariant_vcf_gz.output.vcf,
        vcf_index = rules.prepare_invariant_vcf_gz_index.output.index,
        popmap = rules.generate_pixy_popmap.output.popmap
    output:
        pi = "results/{project}/pixy/{project}.{grouping}.pixy_pi.txt"
    log:
        "logs/{project}/pixy_pi.{grouping}.log"
    benchmark:
        "benchmarks/{project}/pixy_pi.{grouping}.txt"
    params:
        window_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("window_size", 10000),
        output_folder = "results/{project}/pixy/",
        output_prefix = "{project}.{grouping}.pixy"
    conda:
        "../envs/pixy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["runtime"]
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats pi --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder} --output_prefix {params.output_prefix} &> {log}
        """


rule pixy_fst:
    input:
        vcf = rules.prepare_invariant_vcf_gz.output.vcf,
        vcf_index = rules.prepare_invariant_vcf_gz_index.output.index,
        popmap = rules.generate_pixy_popmap.output.popmap
    output:
        fst = "results/{project}/pixy/{project}.{grouping}.pixy_fst.txt"
    log:
        "logs/{project}/pixy_fst.{grouping}.log"
    benchmark:
        "benchmarks/{project}/pixy_fst.{grouping}.txt"
    params:
        window_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("window_size", 10000),
        output_folder = "results/{project}/pixy/",
        output_prefix = "{project}.{grouping}.pixy"
    conda:
        "../envs/pixy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["runtime"]
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats fst --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder} --output_prefix {params.output_prefix} &> {log}
        """


rule pixy_dxy:
    input:
        vcf = rules.prepare_invariant_vcf_gz.output.vcf,
        vcf_index = rules.prepare_invariant_vcf_gz_index.output.index,
        popmap = rules.generate_pixy_popmap.output.popmap
    output:
        dxy = "results/{project}/pixy/{project}.{grouping}.pixy_dxy.txt"
    log:
        "logs/{project}/pixy_dxy.{grouping}.log"
    benchmark:
        "benchmarks/{project}/pixy_dxy.{grouping}.txt"
    params:
        window_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("window_size", 10000),
        output_folder = "results/{project}/pixy/",
        output_prefix = "{project}.{grouping}.pixy"
    conda:
        "../envs/pixy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["runtime"]
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats dxy --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder} --output_prefix {params.output_prefix} &> {log}
        """


rule pixy_pi_summary:
    input:
        pi = rules.pixy_pi.output.pi
    output:
        pi = "results/{project}/pixy/{project}.{grouping}.pixy_pi-summary.txt",
    log:
        "logs/{project}/pixy_pi_summary.{grouping}.log"
    benchmark:
        "benchmarks/{project}/pixy_pi_summary.{grouping}.txt"
    params:
        stat = "pi",
        bootstrap_replicates = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("bootstrap_replicates", 1000)
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy_summary"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy_summary"]["runtime"]
    script:
        "../scripts/pixy_summary.py"


rule pixy_fst_summary:
    input:
        fst = rules.pixy_fst.output.fst
    output:
        fst = "results/{project}/pixy/{project}.{grouping}.pixy_fst-summary.txt"
    log:
        "logs/{project}/pixy_fst_summary.{grouping}.log"
    benchmark:
        "benchmarks/{project}/pixy_fst_summary.{grouping}.txt"
    params:
        stat = "fst",
        bootstrap_replicates = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("bootstrap_replicates", 1000)
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy_summary"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy_summary"]["runtime"]
    script:
        "../scripts/pixy_summary.py"


rule pixy_dxy_summary:
    input:
        dxy = rules.pixy_dxy.output.dxy
    output:
        dxy = "results/{project}/pixy/{project}.{grouping}.pixy_dxy-summary.txt"
    log:
        "logs/{project}/pixy_dxy_summary.{grouping}.log"
    benchmark:
        "benchmarks/{project}/pixy_dxy_summary.{grouping}.txt"
    params:
        stat = "dxy",
        bootstrap_replicates = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("bootstrap_replicates", 1000)
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy_summary"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy_summary"]["runtime"]
    script:
        "../scripts/pixy_summary.py"

# Rule to plot FST heatmap with dendrogram
rule plot_pixy_fst_heatmap:
    input:
        fst_summary = rules.pixy_fst_summary.output.fst
    output:
        pdf = "results/{project}/pixy/plots/{project}.{grouping}.pixy_fst_heatmap.pdf",
        rds = "results/{project}/pixy/plots/{project}.{grouping}.pixy_fst_heatmap.rds"
    log:
        "logs/{project}/plot_pixy_fst_heatmap.{grouping}.log"
    benchmark:
        "benchmarks/{project}/plot_pixy_fst_heatmap.{grouping}.txt"
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
        dxy_summary = rules.pixy_dxy_summary.output.dxy
    output:
        pdf = "results/{project}/pixy/plots/{project}.{grouping}.pixy_dxy_heatmap.pdf",
        rds = "results/{project}/pixy/plots/{project}.{grouping}.pixy_dxy_heatmap.rds"
    log:
        "logs/{project}/plot_pixy_dxy_heatmap.{grouping}.log"
    benchmark:
        "benchmarks/{project}/plot_pixy_dxy_heatmap.{grouping}.txt"
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
        pi_summary = "results/{project}/pixy/{project}.{grouping}.pixy_pi-summary.txt",
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/pixy/plots/{project}.{grouping}.pixy_pi_grouped_by_{color_by}.pdf",
        rds = "results/{project}/pixy/plots/{project}.{grouping}.pixy_pi_grouped_by_{color_by}.rds"
    log:
        "logs/{project}/plot_pixy_pi_barplot_{grouping}_grouped_by_{color_by}.log"
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
        pi_summary = "results/{project}/pixy/{project}.{grouping}.pixy_pi-summary.txt",
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/pixy/plots/{project}.{grouping}.pixy_pi_sorted_by_{color_by}.pdf",
        rds = "results/{project}/pixy/plots/{project}.{grouping}.pixy_pi_sorted_by_{color_by}.rds"
    log:
        "logs/{project}/plot_pixy_pi_barplot_{grouping}_sorted_by_{color_by}.log"
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
        pi_summary = "results/{project}/pixy/{project}.{grouping}.pixy_pi-summary.txt",
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/pixy/plots/{project}.{grouping}.pixy_pi_plain.pdf",
        rds = "results/{project}/pixy/plots/{project}.{grouping}.pixy_pi_plain.rds"
    log:
        "logs/{project}/plot_pixy_pi_barplot_{grouping}_plain.log"
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
        indpopdata = rules.generate_popdata.output.indpopdata,
        # Map plotting requires Site-level populations because
        # plot_pixy_maps.R merges summary population names with Site coords.
        summary = "results/{project}/pixy/{project}.Site.pixy_pi-summary.txt",
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
        point_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("point_size", 3),
        map_outline = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("map_outline", True)
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
