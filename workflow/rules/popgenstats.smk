"""
Population-genetic summary statistics (popgenstats module).

1) Weir & Goudet (2017) dosage Fst/Fis + Ho/Hs on the biallelic SNP VCF.
2) Weir & Cockerham basic.stats, pairwise Fst (+ bootstrap CIs), allelic
   richness, and private alleles on the all-SNP (filtered) VCF.

Groupings reuse generate_pixy_popmap so labels match pixy when columns overlap.
"""

rule popgenstats_stats:
    input:
        vcf = rules.select_biallelic_snps.output.biallelic_vcf,
        popmap = rules.generate_pixy_popmap.output.popmap
    output:
        stats = "results/{project}/popgenstats/{project}.{grouping}.popgenstats_stats.tsv"
    params:
        grouping = lambda wildcards: wildcards.grouping,
        min_individuals = lambda wildcards: _popgenstats_setting(
            wildcards.project, "min_individuals", 2
        )
    log:
        "logs/{project}/popgenstats_stats.{grouping}.log"
    benchmark:
        "benchmarks/{project}/popgenstats_stats.{grouping}.txt"
    conda:
        "../envs/popgenstats.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["popgenstats"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["popgenstats"]["runtime"]
    script:
        "../scripts/popgenstats_stats.R"


rule popgenstats_wc:
    """
    Weir & Cockerham diversity, pairwise Fst + bootstrap CIs, allelic richness,
    and private alleles from the all-SNP filtered VCF.
    """
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf,
        popmap = rules.generate_pixy_popmap.output.popmap
    output:
        overall = "results/{project}/popgenstats/{project}.{grouping}.popgenstats_wc_overall.tsv",
        popstats = "results/{project}/popgenstats/{project}.{grouping}.popgenstats_wc_popstats.tsv",
        pairwise_fst = "results/{project}/popgenstats/{project}.{grouping}.popgenstats_wc_pairwise_fst.tsv",
        pairwise_fst_ci = "results/{project}/popgenstats/{project}.{grouping}.popgenstats_wc_pairwise_fst_ci.tsv"
    params:
        grouping = lambda wildcards: wildcards.grouping,
        min_individuals = lambda wildcards: _popgenstats_setting(
            wildcards.project, "min_individuals", 2
        ),
        nboot = lambda wildcards: _popgenstats_setting(
            wildcards.project, "nboot", 1000
        )
    log:
        "logs/{project}/popgenstats_wc.{grouping}.log"
    benchmark:
        "benchmarks/{project}/popgenstats_wc.{grouping}.txt"
    conda:
        "../envs/popgenstats.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["popgenstats"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["popgenstats"]["runtime"]
    script:
        "../scripts/popgenstats_wc.R"


rule plot_popgenstats_barplot:
    """
    Single-statistic barplot for dosage-based Ho, Hs, Fis, or population-specific Fst.
    """
    input:
        stats = rules.popgenstats_stats.output.stats,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/popgenstats/plots/{project}.{grouping}.popgenstats_{stat}.pdf",
        rds = "results/{project}/popgenstats/plots/{project}.{grouping}.popgenstats_{stat}.rds"
    params:
        stat = lambda wildcards: wildcards.stat,
        grouping = lambda wildcards: wildcards.grouping,
        group_colors = lambda wildcards: _popgenstats_group_setting(
            wildcards.project, wildcards.grouping, "colors"
        ),
        population_sort_by = lambda wildcards: _popgenstats_group_setting(
            wildcards.project, wildcards.grouping, "sort_by"
        ),
        width = lambda wildcards: _popgenstats_plot_style_params(wildcards)["width"],
        height = lambda wildcards: _popgenstats_plot_style_params(wildcards)["height"],
        axis_title_size = lambda wildcards: _popgenstats_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size = lambda wildcards: _popgenstats_plot_style_params(wildcards)["axis_text_size"],
    wildcard_constraints:
        stat="Ho|Hs|Fis|Fst"
    log:
        "logs/{project}/plot_popgenstats_barplot.{grouping}.{stat}.log"
    benchmark:
        "benchmarks/{project}/plot_popgenstats_barplot.{grouping}.{stat}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_popgenstats_barplot.R"


rule plot_popgenstats_wc_barplot:
    """
    Barplot for rarefied allelic richness or private-allele counts (WC popstats).
    """
    input:
        stats = rules.popgenstats_wc.output.popstats,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/popgenstats/plots/{project}.{grouping}.popgenstats_wc_{stat}.pdf",
        rds = "results/{project}/popgenstats/plots/{project}.{grouping}.popgenstats_wc_{stat}.rds"
    params:
        stat = lambda wildcards: wildcards.stat,
        grouping = lambda wildcards: wildcards.grouping,
        group_colors = lambda wildcards: _popgenstats_group_setting(
            wildcards.project, wildcards.grouping, "colors"
        ),
        population_sort_by = lambda wildcards: _popgenstats_group_setting(
            wildcards.project, wildcards.grouping, "sort_by"
        ),
        width = lambda wildcards: _popgenstats_plot_style_params(wildcards)["width"],
        height = lambda wildcards: _popgenstats_plot_style_params(wildcards)["height"],
        axis_title_size = lambda wildcards: _popgenstats_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size = lambda wildcards: _popgenstats_plot_style_params(wildcards)["axis_text_size"],
    wildcard_constraints:
        stat="Ar|private_alleles"
    log:
        "logs/{project}/plot_popgenstats_wc_barplot.{grouping}.{stat}.log"
    benchmark:
        "benchmarks/{project}/plot_popgenstats_wc_barplot.{grouping}.{stat}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_popgenstats_barplot.R"


rule plot_popgenstats_stats:
    """
    Stacked Ho / Hs / Fis / population-specific Fst panels (pixy diversity-panel style).
    """
    input:
        stats = rules.popgenstats_stats.output.stats,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/popgenstats/plots/{project}.{grouping}.popgenstats_stats.pdf",
        rds = "results/{project}/popgenstats/plots/{project}.{grouping}.popgenstats_stats.rds"
    params:
        grouping = lambda wildcards: wildcards.grouping,
        group_colors = lambda wildcards: _popgenstats_group_setting(
            wildcards.project, wildcards.grouping, "colors"
        ),
        population_sort_by = lambda wildcards: _popgenstats_group_setting(
            wildcards.project, wildcards.grouping, "sort_by"
        ),
        width = lambda wildcards: _popgenstats_plot_style_params(wildcards)["width"],
        height = lambda wildcards: _popgenstats_plot_style_params(wildcards)["height"],
        axis_title_size = lambda wildcards: _popgenstats_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size = lambda wildcards: _popgenstats_plot_style_params(wildcards)["axis_text_size"],
    log:
        "logs/{project}/plot_popgenstats_stats.{grouping}.log"
    benchmark:
        "benchmarks/{project}/plot_popgenstats_stats.{grouping}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_popgenstats_stats.R"


rule plot_popgenstats_pairwise_fst_heatmap:
    """
    Heatmap of pairwise Weir & Cockerham Fst with hierarchical clustering.
    """
    input:
        fst_summary = rules.popgenstats_wc.output.pairwise_fst
    output:
        pdf = "results/{project}/popgenstats/plots/{project}.{grouping}.popgenstats_wc_pairwise_fst_heatmap.pdf",
        rds = "results/{project}/popgenstats/plots/{project}.{grouping}.popgenstats_wc_pairwise_fst_heatmap.rds"
    params:
        width = lambda wildcards: _popgenstats_plot_style_params(wildcards)["width"],
        height = lambda wildcards: _popgenstats_plot_style_params(wildcards)["height"],
        axis_title_size = lambda wildcards: _popgenstats_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size = lambda wildcards: _popgenstats_plot_style_params(wildcards)["axis_text_size"],
    log:
        "logs/{project}/plot_popgenstats_pairwise_fst_heatmap.{grouping}.log"
    benchmark:
        "benchmarks/{project}/plot_popgenstats_pairwise_fst_heatmap.{grouping}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_popgenstats_pairwise_fst_heatmap.R"
