"""
Rules for calculating and plotting missingness statistics (individual and locus level).
Missingness is calculated for both:
1. Filtered VCF (before thinning) - {project}.filtered.imiss/lmiss
2. Thinned VCF (after thinning) - {project}.biallelic_snps_thinned.imiss/lmiss
"""

# Missingness from filtered VCF (before thinning)
rule calculate_missing_indv_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        imiss = "results/{project}/missingness_data/filtered/{project}.filtered.imiss"
    log:
        "logs/{project}/calculate_missing_indv_filtered.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/missingness_data/filtered/{wildcards.project}.filtered"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        mkdir -p results/{wildcards.project}/missingness_data/filtered
        vcftools --gzvcf {input.vcf} \
                 --missing-indv \
                 --out {params.out_prefix} &> {log}
        """

rule calculate_missing_loci_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        lmiss = "results/{project}/missingness_data/filtered/{project}.filtered.lmiss"
    log:
        "logs/{project}/calculate_missing_loci_filtered.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/missingness_data/filtered/{wildcards.project}.filtered"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        mkdir -p results/{wildcards.project}/missingness_data/filtered
        vcftools --gzvcf {input.vcf} \
                 --missing-site \
                 --out {params.out_prefix} &> {log}
        mv {params.out_prefix}.lmiss {output.lmiss} || cp {params.out_prefix}.lmiss {output.lmiss}
        """

# Missingness from thinned VCF (after thinning)
rule calculate_missing_indv_thinned:
    input:
        vcf = rules.thin_vcf.output.vcf
    output:
        imiss = "results/{project}/missingness_data/thinned/{project}.biallelic_snps_thinned.imiss"
    log:
        "logs/{project}/calculate_missing_indv_thinned.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/missingness_data/thinned/{wildcards.project}.biallelic_snps_thinned"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        mkdir -p results/{wildcards.project}/missingness_data/thinned
        vcftools --gzvcf {input.vcf} \
                 --missing-indv \
                 --out {params.out_prefix} &> {log}
        """

rule calculate_missing_loci_thinned:
    input:
        vcf = rules.thin_vcf.output.vcf
    output:
        lmiss = "results/{project}/missingness_data/thinned/{project}.biallelic_snps_thinned.lmiss"
    log:
        "logs/{project}/calculate_missing_loci_thinned.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/missingness_data/thinned/{wildcards.project}.biallelic_snps_thinned"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        mkdir -p results/{wildcards.project}/missingness_data/thinned
        vcftools --gzvcf {input.vcf} \
                 --missing-site \
                 --out {params.out_prefix} &> {log}
        mv {params.out_prefix}.lmiss {output.lmiss} || cp {params.out_prefix}.lmiss {output.lmiss}
        """

# Plot imiss histogram - filtered dataset
rule plot_imiss_histogram_filtered:
    input:
        imiss = rules.calculate_missing_indv_filtered.output.imiss
    output:
        pdf = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_histogram.pdf",
        rds = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_histogram.rds"
    log:
        "logs/{project}/plot_imiss_histogram_filtered.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

# Plot imiss histogram - thinned dataset
rule plot_imiss_histogram_thinned:
    input:
        imiss = rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_histogram.pdf",
        rds = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_histogram.rds"
    log:
        "logs/{project}/plot_imiss_histogram_thinned.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

# Plot lmiss histogram - filtered dataset
rule plot_lmiss_histogram_filtered:
    input:
        lmiss = rules.calculate_missing_loci_filtered.output.lmiss
    output:
        pdf = "results/{project}/missingness_data/filtered/plots/{project}.filtered.lmiss_histogram.pdf",
        rds = "results/{project}/missingness_data/filtered/plots/{project}.filtered.lmiss_histogram.rds"
    log:
        "logs/{project}/plot_lmiss_histogram_filtered.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

# Plot lmiss histogram - thinned dataset
rule plot_lmiss_histogram_thinned:
    input:
        lmiss = rules.calculate_missing_loci_thinned.output.lmiss
    output:
        pdf = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.lmiss_histogram.pdf",
        rds = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.lmiss_histogram.rds"
    log:
        "logs/{project}/plot_lmiss_histogram_thinned.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

# Plot imiss barplot grouped by stratification - filtered dataset
rule plot_imiss_barplot_grouped_filtered:
    input:
        imiss = rules.calculate_missing_indv_filtered.output.imiss,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_barplot-grouped-{color_by}.pdf",
        rds = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_barplot-grouped-{color_by}.rds"
    log:
        "logs/{project}/plot_imiss_barplot_{color_by}_grouped_filtered.log"
    wildcard_constraints:
        color_by="(?!sorted-)(?!plain).*"
    params:
        color_by = lambda wildcards: wildcards.color_by if wildcards.color_by != "none" else None,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        plot_type = "grouped"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_imiss_barplot.R"

# Plot imiss barplot sorted by missingness (colored by stratification) - filtered dataset
rule plot_imiss_barplot_sorted_filtered:
    input:
        imiss = rules.calculate_missing_indv_filtered.output.imiss,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_barplot-sorted-{color_by}.pdf",
        rds = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_barplot-sorted-{color_by}.rds"
    log:
        "logs/{project}/plot_imiss_barplot_{color_by}_sorted_filtered.log"
    wildcard_constraints:
        color_by="(?!grouped-)(?!plain).*"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        plot_type = "sorted"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_imiss_barplot.R"

# Plot imiss barplot plain (no coloring) - filtered dataset
rule plot_imiss_barplot_plain_filtered:
    input:
        imiss = rules.calculate_missing_indv_filtered.output.imiss,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_barplot-plain.pdf",
        rds = "results/{project}/missingness_data/filtered/plots/{project}.filtered.imiss_barplot-plain.rds"
    log:
        "logs/{project}/plot_imiss_barplot_plain_filtered.log"
    params:
        color_by = None,
        pca_colors = None,
        plot_type = "plain"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_imiss_barplot.R"

# Plot imiss barplot grouped by stratification - thinned dataset
rule plot_imiss_barplot_grouped_thinned:
    input:
        imiss = rules.calculate_missing_indv_thinned.output.imiss,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_barplot-grouped-{color_by}.pdf",
        rds = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_barplot-grouped-{color_by}.rds"
    log:
        "logs/{project}/plot_imiss_barplot_{color_by}_grouped_thinned.log"
    wildcard_constraints:
        color_by="(?!sorted-)(?!plain).*"
    params:
        color_by = lambda wildcards: wildcards.color_by if wildcards.color_by != "none" else None,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        plot_type = "grouped"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_imiss_barplot.R"

# Plot imiss barplot sorted by missingness (colored by stratification) - thinned dataset
rule plot_imiss_barplot_sorted_thinned:
    input:
        imiss = rules.calculate_missing_indv_thinned.output.imiss,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_barplot-sorted-{color_by}.pdf",
        rds = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_barplot-sorted-{color_by}.rds"
    log:
        "logs/{project}/plot_imiss_barplot_{color_by}_sorted_thinned.log"
    wildcard_constraints:
        color_by="(?!grouped-)(?!plain).*"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        plot_type = "sorted"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_imiss_barplot.R"

# Plot imiss barplot plain (no coloring) - thinned dataset
rule plot_imiss_barplot_plain_thinned:
    input:
        imiss = rules.calculate_missing_indv_thinned.output.imiss,
        popdata = rules.generate_popdata.output.indpopdata
    output:
        pdf = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_barplot-plain.pdf",
        rds = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.imiss_barplot-plain.rds"
    log:
        "logs/{project}/plot_imiss_barplot_plain_thinned.log"
    params:
        color_by = None,
        pca_colors = None,
        plot_type = "plain"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_imiss_barplot.R"

# Plot lmiss barplot plain - filtered dataset
rule plot_lmiss_barplot_plain_filtered:
    input:
        lmiss = rules.calculate_missing_loci_filtered.output.lmiss
    output:
        pdf = "results/{project}/missingness_data/filtered/plots/{project}.filtered.lmiss_barplot-plain.pdf",
        rds = "results/{project}/missingness_data/filtered/plots/{project}.filtered.lmiss_barplot-plain.rds"
    log:
        "logs/{project}/plot_lmiss_barplot_plain_filtered.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_lmiss_barplot.R"

# Plot lmiss barplot plain - thinned dataset
rule plot_lmiss_barplot_plain_thinned:
    input:
        lmiss = rules.calculate_missing_loci_thinned.output.lmiss
    output:
        pdf = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.lmiss_barplot-plain.pdf",
        rds = "results/{project}/missingness_data/thinned/plots/{project}.biallelic_snps_thinned.lmiss_barplot-plain.rds"
    log:
        "logs/{project}/plot_lmiss_barplot_plain_thinned.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_lmiss_barplot.R"

