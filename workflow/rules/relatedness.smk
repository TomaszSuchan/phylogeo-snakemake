# Whole-sample relatedness: KING and PC-Relate on the full filtered dataset.


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
