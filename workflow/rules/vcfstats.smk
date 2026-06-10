# VCF statistics from original VCF (after user sample subset, before relatedness filtering)
rule calculate_stats_vcf_original:
    input:
        vcf = rules.subset_vcf.output.vcf
    output:
        stats = "results/{project}/stats_vcf/original/{project}.subset.vcf_stats.txt"
    log:
        "logs/{project}/calculate_stats_vcf_original.log"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/calculate_vcf_stats.py"

# VCF statistics from filtered VCF (before thinning)
rule calculate_stats_vcf_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        stats = "results/{project}/stats_vcf/filtered/{project}.filtered.vcf_stats.txt"
    log:
        "logs/{project}/calculate_stats_vcf_filtered.log"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/calculate_vcf_stats.py"

# VCF statistics from thinned VCF (after thinning)
rule calculate_stats_vcf_thinned:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        stats = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.vcf_stats.txt"
    log:
        "logs/{project}/calculate_stats_vcf_thinned.log"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/calculate_vcf_stats.py"

# VCF statistics from the un-thinned biallelic SNP VCF (MAC>1, optional MAF)
rule calculate_stats_vcf_biallelic:
    input:
        vcf = rules.select_biallelic_snps.output.biallelic_vcf
    output:
        stats = "results/{project}/stats_vcf/biallelic/{project}.biallelic_snps.vcf_stats.txt"
    log:
        "logs/{project}/calculate_stats_vcf_biallelic.log"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/calculate_vcf_stats.py"

# Sequencing depth statistics from original VCF (after user sample subset, before relatedness filtering)
rule calculate_depth_vcf_original:
    input:
        vcf = rules.subset_vcf.output.vcf
    output:
        individual_depth = "results/{project}/stats_vcf/original/{project}.subset.idepth",
        site_depth = "results/{project}/stats_vcf/original/{project}.subset.ldepth.mean"
    log:
        "logs/{project}/calculate_depth_vcf_original.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/stats_vcf/original/{wildcards.project}.subset"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --depth \
                 --out {params.out_prefix} > {log} 2>&1
        vcftools --gzvcf {input.vcf} \
                 --site-mean-depth \
                 --out {params.out_prefix} >> {log} 2>&1
        """

rule summarize_depth_vcf_original:
    input:
        individual_depth = rules.calculate_depth_vcf_original.output.individual_depth,
        site_depth = rules.calculate_depth_vcf_original.output.site_depth
    output:
        summary = "results/{project}/stats_vcf/original/{project}.subset.depth_summary.txt"
    log:
        "logs/{project}/summarize_depth_vcf_original.log"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/summarize_depth_stats.py"

rule plot_idepth_histogram_original:
    input:
        idepth = rules.calculate_depth_vcf_original.output.individual_depth
    output:
        pdf = "results/{project}/stats_vcf/original/plots/{project}.subset.idepth_histogram.pdf",
        rds = "results/{project}/stats_vcf/original/plots/{project}.subset.idepth_histogram.rds",
        summary = "results/{project}/stats_vcf/original/{project}.subset.idepth_summary.txt"
    log:
        "logs/{project}/plot_idepth_histogram_original.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_depth_histogram.R"

rule plot_ldepth_histogram_original:
    input:
        ldepth = rules.calculate_depth_vcf_original.output.site_depth
    output:
        pdf = "results/{project}/stats_vcf/original/plots/{project}.subset.ldepth_histogram.pdf",
        rds = "results/{project}/stats_vcf/original/plots/{project}.subset.ldepth_histogram.rds",
        summary = "results/{project}/stats_vcf/original/{project}.subset.ldepth_summary.txt"
    log:
        "logs/{project}/plot_ldepth_histogram_original.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_depth_histogram.R"

# Sequencing depth statistics from filtered VCF (before thinning)
rule calculate_depth_vcf_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        individual_depth = "results/{project}/stats_vcf/filtered/{project}.filtered.idepth",
        site_depth = "results/{project}/stats_vcf/filtered/{project}.filtered.ldepth.mean"
    log:
        "logs/{project}/calculate_depth_vcf_filtered.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/stats_vcf/filtered/{wildcards.project}.filtered"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --depth \
                 --out {params.out_prefix} > {log} 2>&1
        vcftools --gzvcf {input.vcf} \
                 --site-mean-depth \
                 --out {params.out_prefix} >> {log} 2>&1
        """

rule summarize_depth_vcf_filtered:
    input:
        individual_depth = rules.calculate_depth_vcf_filtered.output.individual_depth,
        site_depth = rules.calculate_depth_vcf_filtered.output.site_depth
    output:
        summary = "results/{project}/stats_vcf/filtered/{project}.filtered.depth_summary.txt"
    log:
        "logs/{project}/summarize_depth_vcf_filtered.log"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/summarize_depth_stats.py"

rule plot_idepth_histogram_filtered:
    input:
        idepth = rules.calculate_depth_vcf_filtered.output.individual_depth
    output:
        pdf = "results/{project}/stats_vcf/filtered/plots/{project}.filtered.idepth_histogram.pdf",
        rds = "results/{project}/stats_vcf/filtered/plots/{project}.filtered.idepth_histogram.rds",
        summary = "results/{project}/stats_vcf/filtered/{project}.filtered.idepth_summary.txt"
    log:
        "logs/{project}/plot_idepth_histogram_filtered.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_depth_histogram.R"

rule plot_ldepth_histogram_filtered:
    input:
        ldepth = rules.calculate_depth_vcf_filtered.output.site_depth
    output:
        pdf = "results/{project}/stats_vcf/filtered/plots/{project}.filtered.ldepth_histogram.pdf",
        rds = "results/{project}/stats_vcf/filtered/plots/{project}.filtered.ldepth_histogram.rds",
        summary = "results/{project}/stats_vcf/filtered/{project}.filtered.ldepth_summary.txt"
    log:
        "logs/{project}/plot_ldepth_histogram_filtered.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_depth_histogram.R"

# Sequencing depth statistics from thinned VCF (after thinning)
rule calculate_depth_vcf_thinned:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        individual_depth = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.idepth",
        site_depth = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.ldepth.mean"
    log:
        "logs/{project}/calculate_depth_vcf_thinned.log"
    params:
        out_prefix = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --depth \
                 --out {params.out_prefix} > {log} 2>&1
        vcftools --gzvcf {input.vcf} \
                 --site-mean-depth \
                 --out {params.out_prefix} >> {log} 2>&1
        """

rule summarize_depth_vcf_thinned:
    input:
        individual_depth = rules.calculate_depth_vcf_thinned.output.individual_depth,
        site_depth = rules.calculate_depth_vcf_thinned.output.site_depth
    output:
        summary = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.depth_summary.txt"
    log:
        "logs/{project}/summarize_depth_vcf_thinned.log"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/summarize_depth_stats.py"

rule plot_idepth_histogram_thinned:
    input:
        idepth = rules.calculate_depth_vcf_thinned.output.individual_depth
    output:
        pdf = "results/{project}/stats_vcf/thinned/plots/{project}.biallelic_snps.idepth_histogram.pdf",
        rds = "results/{project}/stats_vcf/thinned/plots/{project}.biallelic_snps.idepth_histogram.rds",
        summary = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.idepth_summary.txt"
    log:
        "logs/{project}/plot_idepth_histogram_thinned.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_depth_histogram.R"

rule plot_ldepth_histogram_thinned:
    input:
        ldepth = rules.calculate_depth_vcf_thinned.output.site_depth
    output:
        pdf = "results/{project}/stats_vcf/thinned/plots/{project}.biallelic_snps.ldepth_histogram.pdf",
        rds = "results/{project}/stats_vcf/thinned/plots/{project}.biallelic_snps.ldepth_histogram.rds",
        summary = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.ldepth_summary.txt"
    log:
        "logs/{project}/plot_ldepth_histogram_thinned.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_depth_histogram.R"

# VCF statistics from original VCF (after user sample subset, before relatedness filtering)
rule chromosome_stats_vcf_original:
    input:
        vcf = rules.subset_vcf.output.vcf
    output:
        stats = "results/{project}/stats_vcf/original/{project}.subset.chromosome_stats.txt"
    log:
        "logs/{project}/chromosome_stats_vcf_original.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        zgrep -v "^#" {input.vcf} \
        | cut -f1 \
        | sort \
        | uniq -c \
        | awk '{{print $2, $1}}' > {output.stats} 2> {log}
        """

# VCF statistics from filtered VCF (before thinning)
rule chromosome_stats_vcf_filtered:
    input:
        vcf = rules.subset_vcf_after_relatedness.output.vcf
    output:
        stats = "results/{project}/stats_vcf/filtered/{project}.filtered.chromosome_stats.txt"
    log:
        "logs/{project}/chromosome_stats_vcf_filtered.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        zgrep -v "^#" {input.vcf} \
        | cut -f1 \
        | sort \
        | uniq -c \
        | awk '{{print $2, $1}}' > {output.stats} 2> {log}
        """

# VCF statistics from thinned VCF (after thinning)
rule chromosome_stats_vcf_thinned:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        stats = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.chromosome_stats.txt"
    log:
        "logs/{project}/chromosome_stats_vcf_thinned.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        zgrep -v "^#" {input.vcf} \
        | cut -f1 \
        | sort \
        | uniq -c \
        | awk '{{print $2, $1}}' > {output.stats} 2> {log}
        """

