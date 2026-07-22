"""
GONE2: historical Ne from linkage disequilibrium (https://github.com/esrud/GONE2).
Uses shared per-population VCFs from gone2_currentne2_common prep.
Runs multiple random seeds (-S) and aggregates mean Ne ± SD.
"""


rule gone2_install:
    output:
        sentinel=touch(".snakemake/gone2_installed"),
    params:
        install_dir=".snakemake/gone2",
        root=os.getcwd(),
    log:
        "logs/gone2_install.log"
    conda:
        "../envs/gone2.yaml"
    localrule: True
    shell:
        r"""
        set -euo pipefail
        cd "{params.root}"
        mkdir -p "$(dirname "{log}")" "{params.install_dir}"
        : > "{log}"
        INSTALL_DIR="{params.install_dir}"
        if [ ! -x "$INSTALL_DIR/gone2" ]; then
            rm -rf "$INSTALL_DIR/src"
            git clone --depth 1 https://github.com/esrud/GONE2 "$INSTALL_DIR/src" >> "{log}" 2>&1
            (
                cd "$INSTALL_DIR/src"
                if [ "$(uname -s)" = "Darwin" ]; then
                    make macos || make gone
                else
                    make gone
                fi
                test -x gone2 || {{ echo "ERROR: gone2 binary was not built" >&2; exit 1; }}
                cp -f gone2 "../gone2"
                chmod +x "../gone2"
            ) >> "{log}" 2>&1
        fi
        touch "{output.sentinel}"
        """


rule gone2_run_seed:
    input:
        install=rules.gone2_install.output.sentinel,
        vcf=rules.gone2_currentne2_common_subset_vcf.output.vcf,
        chrom_filter=rules.gone2_currentne2_common_subset_vcf.output.chrom_filter,
    output:
        ne="results/{project}/gone2/seeds/{project}.{stratum}.seed{seed}_GONE2_Ne",
    params:
        gone2_bin=".snakemake/gone2/gone2",
        out_stem=lambda wildcards, output: str(output.ne).removesuffix("_GONE2_Ne"),
        recombination_rate=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("recombination_rate_cM_per_Mb", 2.5),
    log:
        "logs/{project}/gone2_run_seed.{stratum}.seed{seed}.log"
    benchmark:
        "benchmarks/{project}/gone2_run_seed_{stratum}_seed{seed}.txt"
    conda:
        "../envs/gone2.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2", {}).get("threads", config["projects"][wildcards.project]["parameters"]["resources"].get("gone2_currentne2_common", {}).get("threads", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"])),
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2", {}).get("mem_mb", config["projects"][wildcards.project]["parameters"]["resources"].get("gone2_currentne2_common", {}).get("mem_mb", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"])),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2", {}).get("runtime", config["projects"][wildcards.project]["parameters"]["resources"].get("gone2_currentne2_common", {}).get("runtime", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"])),
    shell:
        """
        mkdir -p "$(dirname {output.ne})"
        {params.gone2_bin} -g 0 -r {params.recombination_rate} -t {threads} \
            -S {wildcards.seed} -o {params.out_stem} {input.vcf} > {log} 2>&1
        """


rule gone2_aggregate_seeds:
    input:
        _gone2_seed_ne_inputs,
    output:
        ne="results/{project}/gone2/{project}.{stratum}_GONE2_Ne",
        summary="results/{project}/gone2/{project}.{stratum}_GONE2_Ne_summary.tsv",
        seeds="results/{project}/gone2/{project}.{stratum}_GONE2_Ne_seeds.tsv",
    log:
        "logs/{project}/gone2_aggregate_seeds.{stratum}.log"
    benchmark:
        "benchmarks/{project}/gone2_aggregate_seeds_{stratum}.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/gone2_aggregate_seeds.py"


rule gone2_plot_ne:
    input:
        ne=rules.gone2_aggregate_seeds.output.ne,
        summary=rules.gone2_aggregate_seeds.output.summary,
        seeds=rules.gone2_aggregate_seeds.output.seeds,
    output:
        pdf="results/{project}/gone2/plots/{project}.{stratum}.gone2_ne.pdf",
        rds="results/{project}/gone2/plots/{project}.{stratum}.gone2_ne.rds",
    params:
        width=lambda wildcards: _fig_cm_to_in(config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get("width"), 20.32),
        height=lambda wildcards: _fig_cm_to_in(config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get("height"), 12.7),
        show_seed_trajectories=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get("show_seed_trajectories", True),
    log:
        "logs/{project}/gone2_plot_ne.{stratum}.log"
    benchmark:
        "benchmarks/{project}/gone2_plot_ne_{stratum}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_gone2_ne.R"


rule gone2_plot_ne_combined:
    input:
        unpack(_gone2_ne_inputs),
        populations=rules.gone2_currentne2_common_prepare_samples.output.populations,
    output:
        pdf="results/{project}/gone2/plots/{project}.gone2_ne_combined.pdf",
        rds="results/{project}/gone2/plots/{project}.gone2_ne_combined.rds",
    params:
        out_dir=lambda wildcards: f"results/{wildcards.project}/gone2",
        width=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get(
                "combined_width",
                config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get("width"),
            ),
            25.4,
        ),
        height=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get(
                "combined_height",
                config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get("height"),
            ),
            15.24,
        ),
        legend_title=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get(
            "population_column", "Population"
        ),
        group_colors=lambda wildcards: _gone2_currentne2_common_group_setting(
            wildcards.project,
            config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("population_column", "Site"),
            "colors",
        ),
        group_sort_by=lambda wildcards: _gone2_currentne2_common_group_setting(
            wildcards.project,
            config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("population_column", "Site"),
            "sort_by",
        ),
        show_seed_trajectories=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get(
            "combined_show_seed_trajectories",
            config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get("show_seed_trajectories", True),
        ),
    log:
        "logs/{project}/gone2_plot_ne_combined.log"
    benchmark:
        "benchmarks/{project}/gone2_plot_ne_combined.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_gone2_ne_combined.R"


rule gone2_collect_summary:
    input:
        unpack(_gone2_summary_inputs),
    output:
        summary="results/{project}/gone2/{project}.gone2_summary.tsv",
    params:
        prep_dir=lambda wildcards: f"results/{wildcards.project}/gone2_currentne2_common/vcf",
        out_dir=lambda wildcards: f"results/{wildcards.project}/gone2",
    log:
        "logs/{project}/gone2_collect_summary.log"
    benchmark:
        "benchmarks/{project}/gone2_collect_summary.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/gone2_collect_summary.py"
