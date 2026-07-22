"""
currentNe2: contemporary Ne from linkage disequilibrium (https://github.com/esrud/currentNe2).
Uses shared per-population VCFs from gone2_currentne2_common prep.
"""


rule currentne2_install:
    output:
        sentinel=touch(".snakemake/currentne2_installed"),
    params:
        install_dir=".snakemake/currentne2",
        root=os.getcwd(),
    log:
        "logs/currentne2_install.log"
    conda:
        "../envs/currentne2.yaml"
    # Compile on the submit host; cluster nodes may lack a usable CWD/log layout.
    localrule: True
    shell:
        r"""
        set -euo pipefail
        cd "{params.root}"
        mkdir -p "$(dirname "{log}")" "{params.install_dir}"
        : > "{log}"
        INSTALL_DIR="{params.install_dir}"
        if [ ! -x "$INSTALL_DIR/currentne2" ]; then
            rm -rf "$INSTALL_DIR/src"
            git clone --depth 1 https://github.com/esrud/currentNe2 "$INSTALL_DIR/src" >> "{log}" 2>&1
            (
                cd "$INSTALL_DIR/src"
                make currentne
                test -x currentne2 || {{ echo "ERROR: currentne2 binary was not built" >&2; exit 1; }}
                cp -f currentne2 "../currentne2"
                chmod +x "../currentne2"
            ) >> "{log}" 2>&1
        fi
        touch "{output.sentinel}"
        """


rule currentne2_run:
    input:
        install=rules.currentne2_install.output.sentinel,
        vcf=rules.gone2_currentne2_common_subset_vcf.output.vcf,
        chrom_filter=rules.gone2_currentne2_common_subset_vcf.output.chrom_filter,
    output:
        result="results/{project}/currentne2/{project}.{stratum}_currentNe2_OUTPUT.txt",
    params:
        currentne2_bin=".snakemake/currentne2/currentne2",
        recombination_rate=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("recombination_rate_cM_per_Mb", 2.5),
        extra_flags=lambda wildcards: (
            (f" -k {config['projects'][wildcards.project]['parameters']['currentne2']['k']}"
             if config["projects"][wildcards.project]["parameters"].get("currentne2", {}).get("k", None) is not None
             else "")
            + (f" -s {config['projects'][wildcards.project]['parameters']['currentne2']['n_snps']}"
               if config["projects"][wildcards.project]["parameters"].get("currentne2", {}).get("n_snps", None) is not None
               else "")
            + (" -x" if config["projects"][wildcards.project]["parameters"].get("currentne2", {}).get("metapopulation", False) else "")
        ),
        genome_size_arg=lambda wildcards: (
            f" {config['projects'][wildcards.project]['parameters']['currentne2']['genome_size_M']}"
            if config["projects"][wildcards.project]["parameters"].get("currentne2", {}).get("genome_size_M", None) is not None
            else ""
        ),
    log:
        "logs/{project}/currentne2_run.{stratum}.log"
    benchmark:
        "benchmarks/{project}/currentne2_run_{stratum}.txt"
    conda:
        "../envs/currentne2.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("currentne2", {}).get("threads", config["projects"][wildcards.project]["parameters"]["resources"].get("gone2_currentne2_common", {}).get("threads", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"])),
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("currentne2", {}).get("mem_mb", config["projects"][wildcards.project]["parameters"]["resources"].get("gone2_currentne2_common", {}).get("mem_mb", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"])),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("currentne2", {}).get("runtime", config["projects"][wildcards.project]["parameters"]["resources"].get("gone2_currentne2_common", {}).get("runtime", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"])),
    shell:
        """
        mkdir -p "$(dirname {output.result})"
        {params.currentne2_bin} -r {params.recombination_rate} -t {threads} \
            -o {output.result}{params.extra_flags} {input.vcf}{params.genome_size_arg} > {log} 2>&1
        """


rule currentne2_plot_ne:
    input:
        output=rules.currentne2_run.output.result,
    output:
        pdf="results/{project}/currentne2/plots/{project}.{stratum}.currentne2_ne.pdf",
        rds="results/{project}/currentne2/plots/{project}.{stratum}.currentne2_ne.rds",
    params:
        width=lambda wildcards: _fig_cm_to_in(config["projects"][wildcards.project]["parameters"]["currentne2"].get("plot", {}).get("width"), 20.32),
        height=lambda wildcards: _fig_cm_to_in(config["projects"][wildcards.project]["parameters"]["currentne2"].get("plot", {}).get("height"), 12.7),
    log:
        "logs/{project}/currentne2_plot_ne.{stratum}.log"
    benchmark:
        "benchmarks/{project}/currentne2_plot_ne_{stratum}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_currentne2_ne.R"


rule currentne2_plot_ne_combined:
    input:
        unpack(_currentne2_ne_inputs),
        populations=rules.gone2_currentne2_common_prepare_samples.output.populations,
    output:
        pdf="results/{project}/currentne2/plots/{project}.currentne2_ne_combined.pdf",
        rds="results/{project}/currentne2/plots/{project}.currentne2_ne_combined.rds",
    params:
        out_dir=lambda wildcards: f"results/{wildcards.project}/currentne2",
        width=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"]["currentne2"].get("plot", {}).get(
                "combined_width",
                config["projects"][wildcards.project]["parameters"]["currentne2"].get("plot", {}).get("width"),
            ),
            25.4,
        ),
        height=lambda wildcards: _fig_cm_to_in(
            config["projects"][wildcards.project]["parameters"]["currentne2"].get("plot", {}).get(
                "combined_height",
                config["projects"][wildcards.project]["parameters"]["currentne2"].get("plot", {}).get("height"),
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
    log:
        "logs/{project}/currentne2_plot_ne_combined.log"
    benchmark:
        "benchmarks/{project}/currentne2_plot_ne_combined.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_currentne2_ne_combined.R"


rule currentne2_collect_summary:
    input:
        unpack(_currentne2_summary_inputs),
    output:
        summary="results/{project}/currentne2/{project}.currentne2_summary.tsv",
    params:
        prep_dir=lambda wildcards: f"results/{wildcards.project}/gone2_currentne2_common/vcf",
        out_dir=lambda wildcards: f"results/{wildcards.project}/currentne2",
    log:
        "logs/{project}/currentne2_collect_summary.log"
    benchmark:
        "benchmarks/{project}/currentne2_collect_summary.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/currentne2_collect_summary.py"
