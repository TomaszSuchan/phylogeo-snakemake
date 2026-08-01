"""
Stairway Plot 2: nonparametric Ne(t) from the folded site frequency spectrum.
https://github.com/xiaoming-liu/stairway-plot-v2

Consumes the folded SFS under results/{project}/stairwayplot2/easysfs/{grouping}/.
Each population is run in its own work directory:
results/{project}/stairwayplot2/{grouping}/work/{project}.{stratum}/
"""


rule stairwayplot2_install:
    output:
        stairbuilder=".snakemake/stairwayplot2/stairway_plot_es/Stairbuilder.class",
    params:
        url="https://github.com/xiaoming-liu/stairway-plot-v2/raw/master/stairway_plot_v2.2.zip",
    log:
        "logs/stairwayplot2_install.log"
    conda:
        "../envs/stairwayplot2.yaml"
    localrule: True
    shell:
        r"""
        set -euo pipefail
        TMPZIP="$(mktemp "${{TMPDIR:-/tmp}}/stairwayplot2.XXXXXX.zip")"
        curl -fsSL -o "$TMPZIP" "{params.url}" > {log} 2>&1
        unzip -qoj "$TMPZIP" 'stairway_plot_v2.2/stairway_plot_es/*' \
            -d "$(dirname {output.stairbuilder})" >> {log} 2>&1
        rm -f "$TMPZIP"
        """


rule stairwayplot2_blueprint:
    """Stairway Plot 2 blueprint, written into the population work directory."""
    input:
        sfs_dir=rules.easysfs_run.output.sfs_dir,
        L=rules.easysfs_count_L.output.L,
        samples=rules.easysfs_prepare_samples.output.samples_dir,
    output:
        blueprint="results/{project}/stairwayplot2/{grouping}/work/{project}.{stratum}/{project}.{stratum}.blueprint",
        stats="results/{project}/stairwayplot2/{grouping}/{project}.{stratum}.sfs_stats.tsv",
    params:
        proj_file="results/{project}/stairwayplot2/easysfs/{grouping}/samples/{project}.{stratum}.proj.txt",
        exclude_singletons=lambda wildcards: config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("exclude_singletons", False),
        pct_training=lambda wildcards: config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("pct_training", 0.67),
        ninput=lambda wildcards: config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("ninput", 200),
        nrand=lambda wildcards: config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("nrand", None),
        seed=lambda wildcards: config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("seed", None),
        mu=lambda wildcards: config["projects"][wildcards.project]["parameters"]["stairwayplot2"]["mu"],
        year_per_generation=lambda wildcards: config["projects"][wildcards.project]["parameters"]["stairwayplot2"]["year_per_generation"],
        stairway_plot_dir=os.path.abspath(".snakemake/stairwayplot2/stairway_plot_es"),
    log:
        "logs/{project}/stairwayplot2_blueprint.{grouping}.{stratum}.log"
    benchmark:
        "benchmarks/{project}/stairwayplot2_blueprint_{grouping}_{stratum}.txt"
    conda:
        "../envs/stairwayplot2.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/stairwayplot2_blueprint.py"


rule stairwayplot2_run:
    """Stairbuilder writes a batch script; its Step 1 jobs are independent and run in parallel."""
    input:
        stairbuilder=rules.stairwayplot2_install.output.stairbuilder,
        blueprint=rules.stairwayplot2_blueprint.output.blueprint,
    output:
        summary="results/{project}/stairwayplot2/{grouping}/work/{project}.{stratum}/run/{project}_{stratum}.final.summary",
    params:
        stairway_plot_dir=os.path.abspath(".snakemake/stairwayplot2/stairway_plot_es"),
    log:
        "logs/{project}/stairwayplot2_run.{grouping}.{stratum}.log"
    benchmark:
        "benchmarks/{project}/stairwayplot2_run_{grouping}_{stratum}.txt"
    conda:
        "../envs/stairwayplot2.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("stairwayplot2", {}).get("threads", 8),
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("stairwayplot2", {}).get("mem_mb", 16000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("stairwayplot2", {}).get("runtime", 2880),
    shell:
        r"""
        set -euo pipefail
        {{
        cd "$(dirname {input.blueprint})"
        blueprint="$(basename {input.blueprint})"
        java -cp {params.stairway_plot_dir} Stairbuilder "$blueprint"
        grep    'training_testing' "$blueprint.sh" | xargs -P {threads} -I CMD bash -c CMD
        grep -v 'training_testing' "$blueprint.sh" | bash
        }} > {log} 2>&1
        """


rule stairwayplot2_parse_summary:
    input:
        summary=rules.stairwayplot2_run.output.summary,
    output:
        tsv="results/{project}/stairwayplot2/{grouping}/{project}.{stratum}.ne.tsv",
    params:
        year_per_generation=lambda wildcards: config["projects"][wildcards.project]["parameters"]["stairwayplot2"]["year_per_generation"],
    log:
        "logs/{project}/stairwayplot2_parse_summary.{grouping}.{stratum}.log"
    benchmark:
        "benchmarks/{project}/stairwayplot2_parse_summary_{grouping}_{stratum}.txt"
    conda:
        "../envs/stairwayplot2.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/stairwayplot2_parse_summary.py"


rule stairwayplot2_plot_ne:
    input:
        ne=rules.stairwayplot2_parse_summary.output.tsv,
    output:
        pdf="results/{project}/stairwayplot2/{grouping}/plots/{project}.{stratum}.stairwayplot2_ne.pdf",
        rds="results/{project}/stairwayplot2/{grouping}/plots/{project}.{stratum}.stairwayplot2_ne.rds",
        pdf_linear="results/{project}/stairwayplot2/{grouping}/plots/{project}.{stratum}.stairwayplot2_ne_linear.pdf",
        rds_linear="results/{project}/stairwayplot2/{grouping}/plots/{project}.{stratum}.stairwayplot2_ne_linear.rds",
        pdf_xlinear_ylog="results/{project}/stairwayplot2/{grouping}/plots/{project}.{stratum}.stairwayplot2_ne_xlinear_ylog.pdf",
        rds_xlinear_ylog="results/{project}/stairwayplot2/{grouping}/plots/{project}.{stratum}.stairwayplot2_ne_xlinear_ylog.rds",
    params:
        width=lambda wildcards: _fig_cm_to_in(config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("plot", {}).get("width"), 20.32),
        height=lambda wildcards: _fig_cm_to_in(config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("plot", {}).get("height"), 12.7),
    log:
        "logs/{project}/stairwayplot2_plot_ne.{grouping}.{stratum}.log"
    benchmark:
        "benchmarks/{project}/stairwayplot2_plot_ne_{grouping}_{stratum}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_stairwayplot2_ne.R"


rule stairwayplot2_plot_ne_combined:
    input:
        ne=lambda wildcards: _stairwayplot2_ne_files(wildcards.project, wildcards.grouping),
        populations=rules.easysfs_prepare_samples.output.populations,
    output:
        pdf="results/{project}/stairwayplot2/{grouping}/plots/{project}.stairwayplot2_ne_combined.pdf",
        rds="results/{project}/stairwayplot2/{grouping}/plots/{project}.stairwayplot2_ne_combined.rds",
        pdf_linear="results/{project}/stairwayplot2/{grouping}/plots/{project}.stairwayplot2_ne_combined_linear.pdf",
        rds_linear="results/{project}/stairwayplot2/{grouping}/plots/{project}.stairwayplot2_ne_combined_linear.rds",
        pdf_xlinear_ylog="results/{project}/stairwayplot2/{grouping}/plots/{project}.stairwayplot2_ne_combined_xlinear_ylog.pdf",
        rds_xlinear_ylog="results/{project}/stairwayplot2/{grouping}/plots/{project}.stairwayplot2_ne_combined_xlinear_ylog.rds",
    params:
        out_dir=lambda wildcards: f"results/{wildcards.project}/stairwayplot2/{wildcards.grouping}",
        width=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("plot", {}).get(
                "combined_width",
                config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("plot", {}).get("width"),
            ),
            25.4,
        ),
        height=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("plot", {}).get(
                "combined_height",
                config["projects"][wildcards.project]["parameters"]["stairwayplot2"].get("plot", {}).get("height"),
            ),
            15.24,
        ),
        legend_title=lambda wildcards: wildcards.grouping,
        group_colors=lambda wildcards: _stairwayplot2_group_setting(
            wildcards.project,
            wildcards.grouping,
            "colors",
        ),
        group_sort_by=lambda wildcards: _stairwayplot2_group_setting(
            wildcards.project,
            wildcards.grouping,
            "sort_by",
        ),
    log:
        "logs/{project}/stairwayplot2_plot_ne_combined.{grouping}.log"
    benchmark:
        "benchmarks/{project}/stairwayplot2_plot_ne_combined_{grouping}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_stairwayplot2_ne_combined.R"


rule stairwayplot2_collect_summary:
    """One QC row per population; also the target that pulls in every per-population output."""
    input:
        unpack(_stairwayplot2_summary_inputs),
    output:
        summary="results/{project}/stairwayplot2/{grouping}/{project}.stairwayplot2_summary.tsv",
    log:
        "logs/{project}/stairwayplot2_collect_summary.{grouping}.log"
    benchmark:
        "benchmarks/{project}/stairwayplot2_collect_summary_{grouping}.txt"
    conda:
        "../envs/stairwayplot2.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    shell:
        r"""
        set -euo pipefail
        awk 'FNR == 1 && NR > 1 {{ next }} {{ print }}' {input.stats} > {output.summary} 2> {log}
        """
