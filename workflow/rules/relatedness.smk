# Rule to get among sample relatedness unadjusted Ajk statistic based 
# on the method of Yang et al. (2010), doi:10.1038/ng.608).
# Expectation: Ajk is one for an individual with themselves
# and zero for individuals within a populations.
rule relatedness:
    input:
        vcf=rules.thin_vcf.output.vcf
    output:
        "results/{project}/relatedness/{project}.relatedness"
    log:
        "logs/{project}/relatedness.log"
    benchmark:
        "benchmarks/{project}/relatedness.txt"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} --relatedness --stdout > {output} 2> {log}
        """

# Rule to get among sample relatedness unadjusted Ajk statistic based 
# on the method of Manichaikul et al. (2010), doi:10.1093/bioinformatics/btq559.
# Expectation: 1/2 for monozygous twins, 1/4 for full-sibs or parent-offspring,
# 1/8 for second degree, 1/16 for third degree, and 0 for unrelated.
rule relatedness2:
    input:
        # Use Snakemake's automatic file selection with multiple possible inputs
        vcf=rules.thin_vcf.output.vcf
    output:
        "results/{project}/relatedness/{project}.relatedness2"
    log:
        "logs/{project}/relatedness2.log"
    benchmark:
        "benchmarks/{project}/relatedness2.txt"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} --relatedness2 --stdout > {output} 2> {log}
        """

# Rule to calculate pairwise IBD sharing using PLINK --genome method.
# Calculates PI_HAT (proportion of IBD alleles shared) and other IBD statistics.
# PI_HAT values: ~0.5 for monozygous twins, ~0.25 for full-sibs/parent-offspring,
# ~0.125 for 2nd-degree relatives, ~0.0625 for 3rd-degree relatives, ~0 for unrelated.
rule relatedness_genome:
    input:
        bed=rules.vcf_to_plink.output.bed,
        bim=rules.vcf_to_plink.output.bim,
        fam=rules.vcf_to_plink.output.fam
    output:
        "results/{project}/relatedness/{project}.genome"
    log:
        "logs/{project}/relatedness_genome.log"
    benchmark:
        "benchmarks/{project}/relatedness_genome.txt"
    params:
        bfile_prefix="results/{project}/filtered_data/{project}.biallelic_snps_thinned",
        output_prefix="results/{project}/relatedness/{project}"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        plink --bfile {params.bfile_prefix} \
              --genome \
              --out {params.output_prefix} >> {log} 2>&1 || exit 1
        rm -f {params.output_prefix}.log {params.output_prefix}.nosex || true
        """

# Rule to calculate KING kinship on filtered data (for analysis, not filtering)
# Similar to relatedness_genome but uses plink2 --make-king-table
rule relatedness_king:
    input:
        bed=rules.vcf_to_plink.output.bed,
        bim=rules.vcf_to_plink.output.bim,
        fam=rules.vcf_to_plink.output.fam
    output:
        "results/{project}/relatedness/{project}.king"
    log:
        "logs/{project}/relatedness_king.log"
    benchmark:
        "benchmarks/{project}/relatedness_king.txt"
    params:
        bfile_prefix="results/{project}/filtered_data/{project}.biallelic_snps_thinned",
        output_prefix="results/{project}/relatedness/{project}"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        plink2 --bfile {params.bfile_prefix} \
               --make-king-table \
               --out {params.output_prefix} >> {log} 2>&1 || exit 1
        
        # Copy KING table to output location
        if [ -f {params.output_prefix}.kin0 ]; then
            cp {params.output_prefix}.kin0 {output} || exit 1
        elif [ -f {params.output_prefix}.kin ]; then
            cp {params.output_prefix}.kin {output} || exit 1
        else
            echo "Error: KING table file not found" >> {log}
            exit 1
        fi
        
        # Cleanup PLINK side-effect files
        rm -f {params.output_prefix}.kin0 {params.output_prefix}.kin {params.output_prefix}.log || true
        """

# Rule to plot genome network graph (based on PI_HAT from plink --genome)
rule plot_genome_network:
    input:
        genome=rules.relatedness_genome.output,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        pdf="results/{project}/relatedness/plots/{project}.genome_network-{color_by}.pdf",
        rds="results/{project}/relatedness/plots/{project}.genome_network-{color_by}.rds"
    log:
        "logs/{project}/plot_genome_network-{color_by}.log"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        relatedness_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("relatedness_colors", None)
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_related_network_genome.R"

# Rule to plot KING network graph
rule plot_king_network:
    input:
        king=rules.relatedness_king.output,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        pdf="results/{project}/relatedness/plots/{project}.king_network-{color_by}.pdf",
        rds="results/{project}/relatedness/plots/{project}.king_network-{color_by}.rds"
    log:
        "logs/{project}/plot_king_network-{color_by}.log"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        relatedness_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_plot", {}).get("relatedness_colors", None)
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_related_network_king.R"
