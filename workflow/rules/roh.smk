# Rule to index subset VCF
rule bcftools_roh:
    input:
        vcf=rules.subset_vcf_after_relatedness.output.vcf
    output:
        roh="results/{project}/roh/{project}.bcftools_roh.txt"
    log:
        "logs/{project}/bcftools_roh.log"
    benchmark:
        "benchmarks/{project}/bcftools_roh.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    shell:
        """
        bcftools +fill-tags {input.vcf} -- -t AF | \
        bcftools roh --AF-dflt 0.5 -G 30 -Or -o {output.roh} &> {log}
        """

rule roh_summary:
    input:
        roh=rules.bcftools_roh.output.roh,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        summary="results/{project}/roh/{project}.roh_summary.txt",
        per_ind="results/{project}/roh/{project}.roh_per_ind.txt",
        plots=directory("results/{project}/roh/plots")
    params:
        group_by = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("roh", {}).get("group_by", ["Site"])
    log:
        "logs/{project}/{project}.roh_summary.log"
    benchmark:
        "benchmarks/{project}/roh_summary.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/roh_summary.R" 