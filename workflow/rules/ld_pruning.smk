"""
Rules for LD pruning VCF files using PLINK.
This module breaks down the LD pruning process into separate steps for easier debugging.
"""

# Step 1: Convert VCF to PLINK format
rule ld_prune_convert_to_plink:
    input:
        vcf = rules.select_biallelic_snps.output.biallelic_vcf
    output:
        bed = temporary("results/{project}/filtered_data/{project}.ld_prune_temp.bed"),
        bim = temporary("results/{project}/filtered_data/{project}.ld_prune_temp.bim"),
        fam = temporary("results/{project}/filtered_data/{project}.ld_prune_temp.fam")
    log:
        "logs/{project}/ld_prune_convert_to_plink.log"
    benchmark:
        "benchmarks/{project}/ld_prune_convert_to_plink.txt"
    params:
        plink_prefix = "results/{project}/filtered_data/{project}.ld_prune_temp"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        plink2 --vcf {input.vcf} \
               --make-bed \
               --out {params.plink_prefix} \
               --allow-extra-chr \
               --double-id >> {log} 2>&1
        """

# Step 2: Perform LD pruning (generate prune.in and prune.out files)
rule ld_prune_calculate:
    input:
        bed = rules.ld_prune_convert_to_plink.output.bed,
        bim = rules.ld_prune_convert_to_plink.output.bim,
        fam = rules.ld_prune_convert_to_plink.output.fam
    output:
        prune_in = temporary("results/{project}/filtered_data/{project}.ld_prune_temp.prune.in"),
        prune_out = temporary("results/{project}/filtered_data/{project}.ld_prune_temp.prune.out")
    log:
        "logs/{project}/ld_prune_calculate.log"
    benchmark:
        "benchmarks/{project}/ld_prune_calculate.txt"
    params:
        plink_prefix = "results/{project}/filtered_data/{project}.ld_prune_temp",
        r2 = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("ld_pruning", {}).get("r2", 0.5),
        window_size = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("ld_pruning", {}).get("window_size", 50),
        step_size = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("ld_pruning", {}).get("step_size", 5)
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        plink2 --bfile {params.plink_prefix} \
               --indep-pairwise {params.window_size} {params.step_size} {params.r2} \
               --bad-ld \
               --out {params.plink_prefix} >> {log} 2>&1
        """

# Step 3: Extract pruned SNPs
rule ld_prune_extract:
    input:
        bed = rules.ld_prune_convert_to_plink.output.bed,
        bim = rules.ld_prune_convert_to_plink.output.bim,
        fam = rules.ld_prune_convert_to_plink.output.fam,
        prune_in = rules.ld_prune_calculate.output.prune_in
    output:
        bed = temporary("results/{project}/filtered_data/{project}.ld_prune_temp.pruned.bed"),
        bim = temporary("results/{project}/filtered_data/{project}.ld_prune_temp.pruned.bim"),
        fam = temporary("results/{project}/filtered_data/{project}.ld_prune_temp.pruned.fam")
    log:
        "logs/{project}/ld_prune_extract.log"
    benchmark:
        "benchmarks/{project}/ld_prune_extract.txt"
    params:
        plink_prefix = "results/{project}/filtered_data/{project}.ld_prune_temp",
        pruned_prefix = "results/{project}/filtered_data/{project}.ld_prune_temp.pruned"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        plink2 --bfile {params.plink_prefix} \
               --extract {input.prune_in} \
               --make-bed \
               --out {params.pruned_prefix} >> {log} 2>&1
        """

# Step 4: Convert pruned PLINK data back to VCF
rule ld_prune_convert_to_vcf:
    input:
        bed = rules.ld_prune_extract.output.bed,
        bim = rules.ld_prune_extract.output.bim,
        fam = rules.ld_prune_extract.output.fam
    output:
        vcf = "results/{project}/filtered_data/{project}.biallelic_snps_ld_pruned.vcf"
    log:
        "logs/{project}/ld_prune_convert_to_vcf.log"
    benchmark:
        "benchmarks/{project}/ld_prune_convert_to_vcf.txt"
    params:
        pruned_prefix = "results/{project}/filtered_data/{project}.ld_prune_temp.pruned",
        output_prefix = lambda wildcards: "results/{project}/filtered_data/{project}.biallelic_snps_ld_pruned".format(project=wildcards.project)
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        plink2 --bfile {params.pruned_prefix} \
               --export vcf \
               --out {params.output_prefix} >> {log} 2>&1
        """

# Step 5: Compress VCF to final output
rule ld_prune_compress:
    input:
        vcf = rules.ld_prune_convert_to_vcf.output.vcf
    output:
        vcf = "results/{project}/filtered_data/{project}.biallelic_snps_ld_pruned.vcf.gz"
    log:
        "logs/{project}/ld_prune_compress.log"
    benchmark:
        "benchmarks/{project}/ld_prune_compress.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        bgzip {input.vcf}
        """
