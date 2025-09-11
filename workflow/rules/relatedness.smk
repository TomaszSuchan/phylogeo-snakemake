# Rule to get among sample relatedness unadjusted Ajk statistic based 
# on the method of Yang et al. (2010), doi:10.1038/ng.608).
# Expectation: Ajk is one for an individual with themselves
# and zero for individuals within a populations.
rule relatedness:
    input:
        # Use Snakemake's automatic file selection with multiple possible inputs
        vcf=rules.sort_vcf.output.vcf
    output:
        vcf= config["analysis_name"] + "/relatedness/out.relatedness",
    conda:
        "../envs/vcftools.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} --relatedness
        """

# Rule to get among sample relatedness unadjusted Ajk statistic based 
# on the method of Manichaikul et al. (2010), doi:10.1093/bioinformatics/btq559.
# Expectation: 1/2 for monozygous twins, 1/4 for full-sibs or parent-offspring,
# 1/8 for d degree, 1/16 for third degree, and 0 for unrelated.
rule relatedness2:
    input:
        # Use Snakemake's automatic file selection with multiple possible inputs
        vcf=rules.sort_vcf.output.vcf
    output:
        vcf= config["analysis_name"] + "/relatedness/out.relatedness2",
    conda:
        "../envs/vcftools.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} --relatedness2
        """
