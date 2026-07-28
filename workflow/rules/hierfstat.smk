"""
Population-specific Fst/Fis (Weir & Goudet 2017) and Ho/Hs via hierfstat.

Uses the biallelic-SNP VCF from `select_biallelic_snps` and reuses the pixy popmap
rule so groupings and population labels stay identical between the two modules.
"""

rule hierfstat_stats:
    input:
        vcf = rules.select_biallelic_snps.output.biallelic_vcf,
        popmap = rules.generate_pixy_popmap.output.popmap
    output:
        stats = "results/{project}/hierfstat/{project}.{grouping}.hierfstat_stats.tsv"
    params:
        grouping = lambda wildcards: wildcards.grouping,
        min_individuals = lambda wildcards: _hierfstat_setting(
            wildcards.project, "min_individuals", 2
        )
    log:
        "logs/{project}/hierfstat_stats.{grouping}.log"
    benchmark:
        "benchmarks/{project}/hierfstat_stats.{grouping}.txt"
    conda:
        "../envs/hierfstat.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["hierfstat"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["hierfstat"]["runtime"]
    script:
        "../scripts/hierfstat_stats.R"


rule plot_hierfstat_barplot:
    """
    Single-statistic barplot for Ho, Hs, Fis, or population-specific Fst.
    """
    input:
        stats = rules.hierfstat_stats.output.stats,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/hierfstat/plots/{project}.{grouping}.hierfstat_{stat}.pdf",
        rds = "results/{project}/hierfstat/plots/{project}.{grouping}.hierfstat_{stat}.rds"
    params:
        stat = lambda wildcards: wildcards.stat,
        grouping = lambda wildcards: wildcards.grouping,
        group_colors = lambda wildcards: _hierfstat_group_setting(
            wildcards.project, wildcards.grouping, "colors"
        ),
        population_sort_by = lambda wildcards: _hierfstat_group_setting(
            wildcards.project, wildcards.grouping, "sort_by"
        ),
        width = lambda wildcards: _hierfstat_plot_style_params(wildcards)["width"],
        height = lambda wildcards: _hierfstat_plot_style_params(wildcards)["height"],
        axis_title_size = lambda wildcards: _hierfstat_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size = lambda wildcards: _hierfstat_plot_style_params(wildcards)["axis_text_size"],
    wildcard_constraints:
        stat="Ho|Hs|Fis|Fst"
    log:
        "logs/{project}/plot_hierfstat_barplot.{grouping}.{stat}.log"
    benchmark:
        "benchmarks/{project}/plot_hierfstat_barplot.{grouping}.{stat}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_hierfstat_barplot.R"


rule plot_hierfstat_stats:
    """
    Stacked Ho / Hs / Fis / population-specific Fst panels (pixy diversity-panel style).
    """
    input:
        stats = rules.hierfstat_stats.output.stats,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/hierfstat/plots/{project}.{grouping}.hierfstat_stats.pdf",
        rds = "results/{project}/hierfstat/plots/{project}.{grouping}.hierfstat_stats.rds"
    params:
        grouping = lambda wildcards: wildcards.grouping,
        group_colors = lambda wildcards: _hierfstat_group_setting(
            wildcards.project, wildcards.grouping, "colors"
        ),
        population_sort_by = lambda wildcards: _hierfstat_group_setting(
            wildcards.project, wildcards.grouping, "sort_by"
        ),
        width = lambda wildcards: _hierfstat_plot_style_params(wildcards)["width"],
        height = lambda wildcards: _hierfstat_plot_style_params(wildcards)["height"],
        axis_title_size = lambda wildcards: _hierfstat_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size = lambda wildcards: _hierfstat_plot_style_params(wildcards)["axis_text_size"],
    log:
        "logs/{project}/plot_hierfstat_stats.{grouping}.log"
    benchmark:
        "benchmarks/{project}/plot_hierfstat_stats.{grouping}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_hierfstat_stats.R"
