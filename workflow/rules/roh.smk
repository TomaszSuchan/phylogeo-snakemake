# Rule to index subset VCF
rule bcftools_roh:
    input:
        vcf=rules.subset_vcf_after_relatedness.output.vcf
    output:
        roh="results/{project}/roh/{project}.bcftools_roh.txt"
    log:
        "logs/{project}/bcftools_roh.log"
    benchmark:
        "benchmarks/{project}/bcftools_roh.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    shell:
        """
        set -euo pipefail
        bcftools +fill-tags {input.vcf} -- -t AF | \
        bcftools roh --AF-dflt 0.5 -G 30 -Or -o {output.roh} &> {log}
        """

rule roh_summary:
    input:
        roh=rules.bcftools_roh.output.roh,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        summary="results/{project}/roh/{project}.roh_summary.txt",
        per_ind="results/{project}/roh/{project}.roh_per_ind.txt",
        stats="results/{project}/roh/{project}.roh_stats_comparisons.txt"
    params:
        group_by = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("roh", {}).get("group_by", ["Site"])
    log:
        "logs/{project}/{project}.roh_summary.log"
    benchmark:
        "benchmarks/{project}/roh_summary.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/roh_summary.R"


rule roh_plots:
    input:
        roh=rules.bcftools_roh.output.roh,
        per_ind=rules.roh_summary.output.per_ind
    output:
        froh_histogram="results/{project}/roh/plots/{project}.froh_histogram.pdf",
        froh_histogram_rds="results/{project}/roh/plots/{project}.froh_histogram.rds",
        nroh_histogram="results/{project}/roh/plots/{project}.nroh_histogram.pdf",
        nroh_histogram_rds="results/{project}/roh/plots/{project}.nroh_histogram.rds",
        length_histogram="results/{project}/roh/plots/{project}.length_histogram.pdf",
        length_histogram_rds="results/{project}/roh/plots/{project}.length_histogram.rds",
        length_segments_histogram="results/{project}/roh/plots/{project}.length_segments_histogram.pdf",
        length_segments_histogram_rds="results/{project}/roh/plots/{project}.length_segments_histogram.rds",
        froh_classes="results/{project}/roh/plots/{project}.froh_classes.pdf",
        froh_classes_rds="results/{project}/roh/plots/{project}.froh_classes.rds",
        froh_classes_histogram="results/{project}/roh/plots/{project}.froh_classes_histogram.pdf",
        froh_classes_histogram_rds="results/{project}/roh/plots/{project}.froh_classes_histogram.rds",
        nroh_classes="results/{project}/roh/plots/{project}.nroh_classes.pdf",
        nroh_classes_rds="results/{project}/roh/plots/{project}.nroh_classes.rds",
        nroh_classes_histogram="results/{project}/roh/plots/{project}.nroh_classes_histogram.pdf",
        nroh_classes_histogram_rds="results/{project}/roh/plots/{project}.nroh_classes_histogram.rds",
        froh_vs_nroh="results/{project}/roh/plots/{project}.froh_vs_nroh.pdf",
        froh_vs_nroh_rds="results/{project}/roh/plots/{project}.froh_vs_nroh.rds"
    params:
        group_col=None,
        width=lambda wildcards: _roh_plot_style_params(wildcards)["width"],
        height=lambda wildcards: _roh_plot_style_params(wildcards)["height"],
        axis_title_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_text_size"],
        point_size=lambda wildcards: _roh_plot_style_params(wildcards)["point_size"],
    log:
        "logs/{project}/{project}.roh_plots.log"
    benchmark:
        "benchmarks/{project}/roh_plots.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/roh_plots.R"


rule roh_group_plots:
    input:
        per_ind=rules.roh_summary.output.per_ind
    output:
        froh="results/{project}/roh/plots/{project}.froh_by_{group_col}.pdf",
        froh_rds="results/{project}/roh/plots/{project}.froh_by_{group_col}.rds",
        froh_classes="results/{project}/roh/plots/{project}.froh_by_{group_col}_classes.pdf",
        froh_classes_rds="results/{project}/roh/plots/{project}.froh_by_{group_col}_classes.rds",
        nroh="results/{project}/roh/plots/{project}.nroh_by_{group_col}.pdf",
        nroh_rds="results/{project}/roh/plots/{project}.nroh_by_{group_col}.rds",
        nroh_classes="results/{project}/roh/plots/{project}.nroh_by_{group_col}_classes.pdf",
        nroh_classes_rds="results/{project}/roh/plots/{project}.nroh_by_{group_col}_classes.rds"
    params:
        group_col=lambda wildcards: wildcards.group_col,
        group_colors=lambda wildcards: _roh_group_setting(
            wildcards.project, wildcards.group_col, "colors"
        ),
        group_sort_by=lambda wildcards: _roh_group_setting(
            wildcards.project, wildcards.group_col, "sort_by"
        ),
        width=lambda wildcards: _roh_plot_style_params(wildcards)["width"],
        height=lambda wildcards: _roh_plot_style_params(wildcards)["height"],
        axis_title_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_text_size"],
        point_size=lambda wildcards: _roh_plot_style_params(wildcards)["point_size"],
    log:
        "logs/{project}/roh_group_plots.{group_col}.log"
    benchmark:
        "benchmarks/{project}/roh_group_plots.{group_col}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/roh_plots.R"


rule plot_froh_map:
    input:
        unpack(_froh_map_inputs),
    output:
        pdf = "results/{project}/roh/plots/{project}.froh_map.pdf",
        rds = "results/{project}/roh/plots/{project}.froh_map.rds"
    params:
        lambda wildcards: _froh_map_params(wildcards)
    log:
        "logs/{project}/plot_froh_map.log"
    benchmark:
        "benchmarks/{project}/plot_froh_map.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/plot_froh_map.R"
