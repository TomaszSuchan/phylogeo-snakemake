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
        stats="results/{project}/roh/{project}.roh_stats_comparisons.txt",
        froh_histogram="results/{project}/roh/plots/{project}.froh_histogram.pdf",
        n_roh_segments_histogram="results/{project}/roh/plots/{project}.n_roh_segments_histogram.pdf",
        total_roh_length_histogram="results/{project}/roh/plots/{project}.total_roh_length_histogram.pdf",
        roh_segment_length_distribution="results/{project}/roh/plots/{project}.roh_segment_length_distribution.pdf",
        roh_by_class="results/{project}/roh/plots/{project}.roh_by_class.pdf",
        froh_by_class_panel="results/{project}/roh/plots/{project}.froh_per_class_panel.pdf",
        froh_vs_n_segments="results/{project}/roh/plots/{project}.froh_vs_n_segments.pdf"
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


rule roh_group_plots:
    input:
        roh=rules.bcftools_roh.output.roh,
        per_ind=rules.roh_summary.output.per_ind,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        froh="results/{project}/roh/plots/{project}.froh_by_{group_col}.pdf",
        froh_class_panel="results/{project}/roh/plots/{project}.froh_by_{group_col}_per_class_panel.pdf",
        nseg="results/{project}/roh/plots/{project}.n_roh_segments_by_{group_col}.pdf",
        roh_class="results/{project}/roh/plots/{project}.roh_class_by_{group_col}.pdf"
    params:
        group_col = lambda wildcards: wildcards.group_col
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
        "../scripts/roh_group_plots.R"


rule plot_roh_froh_map:
    input:
        unpack(_roh_froh_map_inputs),
    output:
        pdf = "results/{project}/roh/plots/{project}.roh_froh_map.pdf",
        rds = "results/{project}/roh/plots/{project}.roh_froh_map.rds"
    params:
        lambda wildcards: _roh_froh_map_params(wildcards)
    log:
        "logs/{project}/plot_roh_froh_map.log"
    benchmark:
        "benchmarks/{project}/plot_roh_froh_map.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/plot_roh_froh_map.R"
