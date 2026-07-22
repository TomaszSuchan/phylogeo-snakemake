# Population-stratified relatedness: per-stratum VCF subsets and within-population estimators.


rule install_related:
    output:
        touch(".snakemake/related_installed"),
    conda:
        "../envs/related.yaml"
    log:
        "logs/install_related.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        Rscript --vanilla -e 'lib <- .libPaths()[1]; unlink(list.files(lib, pattern="^00LOCK-", full.names=TRUE), recursive=TRUE, force=TRUE); if (!requireNamespace("related", quietly=TRUE, lib.loc=lib)) remotes::install_url("https://github.com/timothyfrasier/related/raw/master/related_1.0.tar.gz", upgrade="never", dependencies=TRUE, build_vignettes=FALSE, lib=lib); if (!requireNamespace("related", quietly=TRUE, lib.loc=lib)) stop("Failed to install related package")' >> {log} 2>&1
        """


rule relatedness_prepare_samples:
    input:
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        samples_dir=directory("results/{project}/relatedness_by_pop/samples"),
        populations="results/{project}/relatedness_by_pop/{project}.populations.tsv",
    params:
        population_column=lambda wildcards: _relatedness_stratum_column(
            config["projects"][wildcards.project]["parameters"]
        ),
    log:
        "logs/{project}/relatedness_prepare_samples.log"
    benchmark:
        "benchmarks/{project}/relatedness_prepare_samples.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/relatedness_prepare_samples.py"


rule relatedness_subset_vcf:
    input:
        vcf=lambda wildcards: get_filtered_vcf_output(wildcards),
        samples=rules.relatedness_prepare_samples.output.samples_dir,
    output:
        vcf="results/{project}/relatedness_by_pop/vcf/{project}.{stratum}.vcf.gz",
        index="results/{project}/relatedness_by_pop/vcf/{project}.{stratum}.vcf.gz.csi",
    params:
        samples_file=lambda wildcards: (
            f"results/{wildcards.project}/relatedness_by_pop/samples/"
            f"{wildcards.project}.{wildcards.stratum}.samples.txt"
        ),
    log:
        "logs/{project}/relatedness_subset_vcf.{stratum}.log"
    benchmark:
        "benchmarks/{project}/relatedness_subset_vcf_{stratum}.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["related"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["related"]["runtime"],
    shell:
        r"""
        set -euo pipefail
        bcftools view -S {params.samples_file} -Oz -o {output.vcf} {input.vcf} &> {log}
        bcftools index --csi {output.vcf} &>> {log}
        """


# Yang et al. (2010) Ajk genomic relatedness within each population stratum.
rule relatedness:
    input:
        vcf=rules.relatedness_subset_vcf.output.vcf,
    output:
        "results/{project}/relatedness_by_pop/{project}.{stratum}.relatedness",
    log:
        "logs/{project}/relatedness.{stratum}.log"
    benchmark:
        "benchmarks/{project}/relatedness_{stratum}.txt"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    shell:
        """
        vcftools --gzvcf {input.vcf} --relatedness --stdout > {output} 2> {log}
        """


rule relatedness_summary:
    input:
        unpack(lambda wildcards: _relatedness_stratum_summary_inputs(wildcards, "relatedness")),
    output:
        summary="results/{project}/relatedness_by_pop/{project}.relatedness_summary.tsv",
    params:
        kind="relatedness",
    log:
        "logs/{project}/relatedness_summary.log"
    benchmark:
        "benchmarks/{project}/relatedness_summary.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/relatedness_collect_results.py"


# Manichaikul et al. (2010) kinship within each population stratum.
rule relatedness2:
    input:
        vcf=rules.relatedness_subset_vcf.output.vcf,
    output:
        "results/{project}/relatedness_by_pop/{project}.{stratum}.relatedness2",
    log:
        "logs/{project}/relatedness2.{stratum}.log"
    benchmark:
        "benchmarks/{project}/relatedness2_{stratum}.txt"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    shell:
        """
        vcftools --gzvcf {input.vcf} --relatedness2 --stdout > {output} 2> {log}
        """


rule relatedness2_summary:
    input:
        unpack(lambda wildcards: _relatedness_stratum_summary_inputs(wildcards, "relatedness2")),
    output:
        summary="results/{project}/relatedness_by_pop/{project}.relatedness2_summary.tsv",
    params:
        kind="relatedness2",
    log:
        "logs/{project}/relatedness2_summary.log"
    benchmark:
        "benchmarks/{project}/relatedness2_summary.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/relatedness_collect_results.py"


rule relatedness_stratum_plink:
    input:
        vcf=rules.relatedness_subset_vcf.output.vcf,
    output:
        bed="results/{project}/relatedness_by_pop/plink/{project}.{stratum}.bed",
        bim="results/{project}/relatedness_by_pop/plink/{project}.{stratum}.bim",
        fam="results/{project}/relatedness_by_pop/plink/{project}.{stratum}.fam",
    params:
        output_prefix="results/{project}/relatedness_by_pop/plink/{project}.{stratum}",
    log:
        "logs/{project}/relatedness_stratum_plink.{stratum}.log"
    benchmark:
        "benchmarks/{project}/relatedness_stratum_plink_{stratum}.txt"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --out {params.output_prefix} \
              --allow-extra-chr 0 \
              --double-id >> {log} 2>&1 || exit 1
        """


# PLINK --genome IBD sharing within each population stratum.
rule relatedness_genome:
    input:
        bed=rules.relatedness_stratum_plink.output.bed,
        bim=rules.relatedness_stratum_plink.output.bim,
        fam=rules.relatedness_stratum_plink.output.fam,
    output:
        "results/{project}/relatedness_by_pop/{project}.{stratum}.genome",
    log:
        "logs/{project}/relatedness_genome.{stratum}.log"
    benchmark:
        "benchmarks/{project}/relatedness_genome_{stratum}.txt"
    params:
        bfile_prefix="results/{project}/relatedness_by_pop/plink/{project}.{stratum}",
        output_prefix="results/{project}/relatedness_by_pop/{project}.{stratum}",
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    shell:
        """
        plink --bfile {params.bfile_prefix} \
              --genome \
              --out {params.output_prefix} >> {log} 2>&1 || exit 1
        rm -f {params.output_prefix}.log {params.output_prefix}.nosex || true
        """


rule relatedness_genome_summary:
    input:
        unpack(lambda wildcards: _relatedness_stratum_summary_inputs(wildcards, "genome")),
    output:
        summary="results/{project}/relatedness_by_pop/{project}.genome_summary.tsv",
    params:
        kind="genome",
    log:
        "logs/{project}/relatedness_genome_summary.log"
    benchmark:
        "benchmarks/{project}/relatedness_genome_summary.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/relatedness_collect_results.py"


# R related package (Coancestry; Queller & Goodnight Rxy) within each population stratum.
rule related_coancestry:
    input:
        install=rules.install_related.output,
        vcf=rules.relatedness_subset_vcf.output.vcf,
    output:
        relatedness="results/{project}/relatedness_by_pop/{project}.{stratum}.related_coancestry.tsv",
    params:
        population_label=lambda wildcards: _relatedness_population_label(
            wildcards.project, wildcards.stratum
        ),
        estimators=lambda wildcards: config["projects"][wildcards.project]["parameters"].get(
            "related", {}
        ).get("estimators", {"quellergt": 1}),
    log:
        "logs/{project}/related_coancestry.{stratum}.log",
    benchmark:
        "benchmarks/{project}/related_coancestry_{stratum}.txt",
    conda:
        "../envs/related.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["related"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["related"]["runtime"],
    script:
        "../scripts/run_related_coancestry.R"


rule related_coancestry_summary:
    input:
        unpack(lambda wildcards: _relatedness_stratum_summary_inputs(wildcards, "related_coancestry")),
    output:
        summary="results/{project}/relatedness_by_pop/{project}.related_coancestry_summary.tsv",
    params:
        kind="related_coancestry",
    log:
        "logs/{project}/related_coancestry_summary.log",
    benchmark:
        "benchmarks/{project}/related_coancestry_summary.txt",
    conda:
        "../envs/python.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/relatedness_collect_results.py"


rule plot_genome_network:
    input:
        pairwise=rules.relatedness_genome_summary.output.summary,
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        pdf="results/{project}/relatedness_by_pop/plots/{project}.genome_network-{color_by}.pdf",
        rds="results/{project}/relatedness_by_pop/plots/{project}.genome_network-{color_by}.rds",
    log:
        "logs/{project}/plot_genome_network-{color_by}.log",
    params:
        color_by=lambda wildcards: wildcards.color_by,
        relatedness_colors=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("relatedness_colors", None),
        plot_all=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("plot_all", False),
        threshold_profile="pihat",
        weight_column="PI_HAT",
    conda:
        "../envs/r-plot.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_related_network.R"


rule plot_relatedness_network:
    input:
        pairwise=rules.relatedness_summary.output.summary,
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        pdf="results/{project}/relatedness_by_pop/plots/{project}.relatedness_network-{color_by}.pdf",
        rds="results/{project}/relatedness_by_pop/plots/{project}.relatedness_network-{color_by}.rds",
    log:
        "logs/{project}/plot_relatedness_network-{color_by}.log",
    params:
        color_by=lambda wildcards: wildcards.color_by,
        relatedness_colors=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("relatedness_colors", None),
        plot_all=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("plot_all", False),
        threshold_profile="ajk",
        weight_column="RELATEDNESS_AJK",
    conda:
        "../envs/r-plot.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_related_network.R"


rule plot_relatedness2_network:
    input:
        pairwise=rules.relatedness2_summary.output.summary,
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        pdf="results/{project}/relatedness_by_pop/plots/{project}.relatedness2_network-{color_by}.pdf",
        rds="results/{project}/relatedness_by_pop/plots/{project}.relatedness2_network-{color_by}.rds",
    log:
        "logs/{project}/plot_relatedness2_network-{color_by}.log",
    params:
        color_by=lambda wildcards: wildcards.color_by,
        relatedness_colors=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("relatedness_colors", None),
        plot_all=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("plot_all", False),
        threshold_profile="manichaikul",
        weight_column="RELATEDNESS_PHI",
    conda:
        "../envs/r-plot.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_related_network.R"


rule plot_coancestry_network:
    input:
        pairwise=rules.related_coancestry_summary.output.summary,
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        pdf="results/{project}/relatedness_by_pop/plots/{project}.coancestry_network-{color_by}.pdf",
        rds="results/{project}/relatedness_by_pop/plots/{project}.coancestry_network-{color_by}.rds",
    log:
        "logs/{project}/plot_coancestry_network-{color_by}.log",
    params:
        color_by=lambda wildcards: wildcards.color_by,
        relatedness_colors=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("relatedness_colors", None),
        plot_all=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("plot_all", False),
        threshold_profile="manichaikul",
        weight_column="",
    conda:
        "../envs/r-plot.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_related_network.R"
