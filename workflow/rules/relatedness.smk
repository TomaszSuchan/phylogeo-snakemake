# Stratified relatedness infrastructure (shared by all population-stratified estimators).
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
        samples_dir=directory("results/{project}/relatedness/samples"),
        populations="results/{project}/relatedness/{project}.populations.tsv",
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
        vcf="results/{project}/relatedness/vcf/{project}.{stratum}.vcf.gz",
        index="results/{project}/relatedness/vcf/{project}.{stratum}.vcf.gz.csi",
    params:
        samples_file=lambda wildcards: (
            f"results/{wildcards.project}/relatedness/samples/"
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
        "results/{project}/relatedness/{project}.{stratum}.relatedness",
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
        summary="results/{project}/relatedness/{project}.relatedness_summary.tsv",
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
        "results/{project}/relatedness/{project}.{stratum}.relatedness2",
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
        summary="results/{project}/relatedness/{project}.relatedness2_summary.tsv",
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
        bed="results/{project}/relatedness/plink/{project}.{stratum}.bed",
        bim="results/{project}/relatedness/plink/{project}.{stratum}.bim",
        fam="results/{project}/relatedness/plink/{project}.{stratum}.fam",
    params:
        output_prefix="results/{project}/relatedness/plink/{project}.{stratum}",
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
        "results/{project}/relatedness/{project}.{stratum}.genome",
    log:
        "logs/{project}/relatedness_genome.{stratum}.log"
    benchmark:
        "benchmarks/{project}/relatedness_genome_{stratum}.txt"
    params:
        bfile_prefix="results/{project}/relatedness/plink/{project}.{stratum}",
        output_prefix="results/{project}/relatedness/{project}.{stratum}",
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
        summary="results/{project}/relatedness/{project}.genome_summary.tsv",
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


# KING kinship on the full sample (robust to population structure without stratification).
rule relatedness_king:
    input:
        bed=rules.vcf_to_plink.output.bed,
        bim=rules.vcf_to_plink.output.bim,
        fam=rules.vcf_to_plink.output.fam,
    output:
        "results/{project}/relatedness/{project}.king",
    log:
        "logs/{project}/relatedness_king.log"
    benchmark:
        "benchmarks/{project}/relatedness_king.txt"
    params:
        bfile_prefix="results/{project}/filtered_data/{project}.biallelic_snps_thinned",
        output_prefix="results/{project}/relatedness/{project}",
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    shell:
        """
        plink2 --bfile {params.bfile_prefix} \
               --make-king-table \
               --out {params.output_prefix} >> {log} 2>&1 || exit 1

        if [ -f {params.output_prefix}.kin0 ]; then
            cp {params.output_prefix}.kin0 {output} || exit 1
        elif [ -f {params.output_prefix}.kin ]; then
            cp {params.output_prefix}.kin {output} || exit 1
        else
            echo "Error: KING table file not found" >> {log}
            exit 1
        fi

        rm -f {params.output_prefix}.kin0 {params.output_prefix}.kin {params.output_prefix}.log || true
        """


# R related package (Coancestry; Queller & Goodnight Rxy) within each population stratum.
rule related_coancestry:
    input:
        install=rules.install_related.output,
        vcf=rules.relatedness_subset_vcf.output.vcf,
    output:
        relatedness="results/{project}/relatedness/{project}.{stratum}.related_coancestry.tsv",
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
        summary="results/{project}/relatedness/{project}.related_coancestry_summary.tsv",
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


# PC-Relate (GENESIS): PC-AiR + PC-Relate on the full sample; ancestry PCs adjust for structure.
rule pcrelate_analysis:
    input:
        vcf=lambda wildcards: get_filtered_vcf_output(wildcards),
    output:
        kinship="results/{project}/relatedness/{project}.pcrelate_kinship.tsv",
    params:
        n_pcs=lambda wildcards: config["projects"][wildcards.project]["parameters"].get(
            "pcrelate", {}
        ).get("n_pcs", 2),
        ld_r2=lambda wildcards: config["projects"][wildcards.project]["parameters"].get(
            "pcrelate", {}
        ).get("ld_r2", 0.2),
        ld_window=lambda wildcards: config["projects"][wildcards.project]["parameters"].get(
            "pcrelate", {}
        ).get("ld_window", 500),
        maf=lambda wildcards: config["projects"][wildcards.project]["parameters"].get(
            "pcrelate", {}
        ).get("maf", 0.01),
        return_ibd_probs=lambda wildcards: config["projects"][wildcards.project]["parameters"].get(
            "pcrelate", {}
        ).get("return_ibd_probs", True),
    log:
        "logs/{project}/pcrelate_analysis.log",
    benchmark:
        "benchmarks/{project}/pcrelate_analysis.txt",
    conda:
        "../envs/pcrelate.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcrelate"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcrelate"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcrelate"]["runtime"],
    script:
        "../scripts/run_pcrelate.R"


rule plot_genome_network:
    input:
        pairwise=rules.relatedness_genome_summary.output.summary,
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        pdf="results/{project}/relatedness/plots/{project}.genome_network-{color_by}.pdf",
        rds="results/{project}/relatedness/plots/{project}.genome_network-{color_by}.rds",
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


rule plot_king_network:
    input:
        pairwise=rules.relatedness_king.output,
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        pdf="results/{project}/relatedness/plots/{project}.king_network-{color_by}.pdf",
        rds="results/{project}/relatedness/plots/{project}.king_network-{color_by}.rds",
    log:
        "logs/{project}/plot_king_network-{color_by}.log",
    params:
        color_by=lambda wildcards: wildcards.color_by,
        relatedness_colors=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("relatedness_colors", None),
        plot_all=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("plot_all", False),
        threshold_profile="king",
        weight_column="KINSHIP",
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
        pdf="results/{project}/relatedness/plots/{project}.relatedness_network-{color_by}.pdf",
        rds="results/{project}/relatedness/plots/{project}.relatedness_network-{color_by}.rds",
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
        pdf="results/{project}/relatedness/plots/{project}.relatedness2_network-{color_by}.pdf",
        rds="results/{project}/relatedness/plots/{project}.relatedness2_network-{color_by}.rds",
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
        pdf="results/{project}/relatedness/plots/{project}.coancestry_network-{color_by}.pdf",
        rds="results/{project}/relatedness/plots/{project}.coancestry_network-{color_by}.rds",
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


rule plot_pcrelate_network:
    input:
        pairwise=rules.pcrelate_analysis.output.kinship,
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        pdf="results/{project}/relatedness/plots/{project}.pcrelate_network-{color_by}.pdf",
        rds="results/{project}/relatedness/plots/{project}.pcrelate_network-{color_by}.rds",
    log:
        "logs/{project}/plot_pcrelate_network-{color_by}.log",
    params:
        color_by=lambda wildcards: wildcards.color_by,
        relatedness_colors=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("relatedness_colors", None),
        plot_all=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("plot_all", False),
        threshold_profile="king",
        weight_column="kin",
    conda:
        "../envs/r-plot.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_related_network.R"
