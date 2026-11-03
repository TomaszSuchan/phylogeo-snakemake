"""
rapidnj-snp: neighbour-joining from variable sites only (filtered VCF → FASTA → RapidNJ).

Runs in parallel with loci-based njtree when both analyses are enabled.
"""


def _rapidnj_snp_params(wildcards):
    params = config["projects"][wildcards.project]["parameters"]
    nj = params.get("njtree", {})
    snp = params.get("rapidnj_snp", {})
    return {
        "distance_model": snp.get("distance_model", nj.get("distance_model", "kim")),
        "bootstraps": int(snp.get("bootstraps", nj.get("bootstraps", 1000))),
        "support_threshold": snp.get("support_threshold", nj.get("support_threshold", 70)),
    }


def _rapidnj_snp_resources(wildcards):
    return config["projects"][wildcards.project]["parameters"]["resources"].get("rapidnj_snp", {})


rule vcf_to_snp_alignment:
    """Encode biallelic SNP genotypes from the analysis VCF as FASTA/PHYLIP."""
    input:
        vcf=lambda wildcards: get_filtered_vcf_output(wildcards),
    output:
        fasta="results/{project}/rapidnj-snp/{project}.snp.fasta",
        phy="results/{project}/rapidnj-snp/{project}.snp.phy",
        summary="results/{project}/rapidnj-snp/{project}.snp_alignment_summary.tsv",
    log:
        "logs/{project}/vcf_to_snp_alignment.log",
    benchmark:
        "benchmarks/{project}/vcf_to_snp_alignment.txt",
    conda:
        "../envs/python.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: _rapidnj_snp_resources(wildcards).get("mem_mb", 16000),
        runtime=lambda wildcards: _rapidnj_snp_resources(wildcards).get("convert_runtime", 120),
    script:
        "../scripts/vcf_to_snp_alignment.py"


rule rapidnj_snp:
    """Neighbour-joining tree with RapidNJ from SNP-only FASTA alignment."""
    input:
        fasta=rules.vcf_to_snp_alignment.output.fasta,
    output:
        treefile="results/{project}/rapidnj-snp/{project}.rapidnj-snp.treefile",
    log:
        "logs/{project}/rapidnj_snp.log",
    benchmark:
        "benchmarks/{project}/rapidnj_snp.txt",
    params:
        distance_model=lambda wildcards: _rapidnj_snp_params(wildcards)["distance_model"],
        bootstrap_flag=lambda wildcards: (
            f"-b {_rapidnj_snp_params(wildcards)['bootstraps']}"
            if _rapidnj_snp_params(wildcards)["bootstraps"] > 0
            else ""
        ),
    conda:
        "../envs/rapidnj.yaml",
    threads: lambda wildcards: _rapidnj_snp_resources(wildcards).get("threads", 12),
    resources:
        mem_mb=lambda wildcards: _rapidnj_snp_resources(wildcards).get("mem_mb", 16000),
        runtime=lambda wildcards: _rapidnj_snp_resources(wildcards).get("runtime", 720),
    shell:
        """
        # -a is the sequence distance model (jc/kim). Do not use -d here: that is HDD cache dir.
        rapidnj {input.fasta} -i fa -a {params.distance_model} -c {threads} \
            -x {output.treefile} {params.bootstrap_flag} > {log} 2>&1
        """


rule plot_rapidnj_snp:
    input:
        treefile=rules.rapidnj_snp.output.treefile,
    output:
        unrooted_pdf="results/{project}/rapidnj-snp/plots/{project}.tree_plot.unrooted.pdf",
        unrooted_rds="results/{project}/rapidnj-snp/plots/{project}.tree_plot.unrooted.rds",
        unrooted_tree="results/{project}/rapidnj-snp/{project}.rapidnj-snp.unrooted.treefile",
    benchmark:
        "benchmarks/{project}/plot_rapidnj_snp.txt",
    params:
        support_threshold=lambda wildcards: _rapidnj_snp_params(wildcards)["support_threshold"],
    conda:
        "../envs/r-ggtree.yaml",
    threads: 1,
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("mem_mb", 8000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("default", {}).get("runtime", 10),
    script:
        "../scripts/plot_nj_tree.R"
