rule generate_popmap:
    input:
        vcf=rules.sort_vcf.output.vcf
    output:
        popmap="results/{project}/popmap.txt"
    params:
        popmap=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popmap", ""),
        separator=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popseparator", "-")
    log:
        "logs/{project}/generate_popmap.log"
    benchmark:
        "benchmarks/{project}/generate_popmap.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        r"""
        if [ -n "{params.popmap}" ] && [ -f "{params.popmap}" ]; then
            cp {params.popmap} {output.popmap}
        else
            bcftools query -l {input.vcf} | \
            awk -v sep="{params.separator}" '{{ 
                split($0,a,sep); 
                pop = (length(a)>1) ? a[1] : "UNKNOWN"; 
                print $0 "\t" pop 
            }}' > {output.popmap}
        fi
        """

rule generate_popdata:
    input:
        popmap=rules.generate_popmap.output.popmap
    output:
        indpopdata="results/{project}/indpopdata.txt"
    params:
        popdata=lambda wildcards: config["projects"][wildcards.project]["parameters"]["popdata"],
    log:
        "logs/{project}/generate_popdata.log"
    benchmark:
        "benchmarks/{project}/generate_popdata.txt"
    conda:
        "../envs/pandas.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/generate_popdata.py"



rule calculate_missing_indv:
    input:
        vcf = rules.thin_vcf.output.vcf
    output:
        imiss = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    log:
        "logs/{project}/calculate_missing_indv.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_thinned"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --missing-indv \
                 --out {params.out_prefix} &> {log}
        """

rule calculate_missing_indv_miss:
    input:
        vcf = rules.filter_missing_vcf.output.vcf
    output:
        imiss = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.imiss"
    log:
        "logs/{project}/calculate_missing_indv_miss_{miss}.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_thinned_miss{wildcards.miss}"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --missing-indv \
                 --out {params.out_prefix} &> {log}
        """