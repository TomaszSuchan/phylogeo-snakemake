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
        
        rm {params.output_prefix}.log {params.output_prefix}.nosex
        """
