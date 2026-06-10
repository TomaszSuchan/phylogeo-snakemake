"""
Rule to generate a manuscript-ready methods.md for each project.

Reads the project config (analyses enabled/disabled, all parameters),
data/mainparams (STRUCTURE BURNIN/NUMREPS), the thinned-VCF stats file
(SNP and RAD-loci counts), the thinned-VCF depth summary, and
workflow/envs/*.yaml (software versions).

Only sections for enabled analyses are emitted. All parameter values and
software versions are filled in automatically. Best-K values are not
included — they are results and belong in the Results section.

Output: results/{project}/methods.md
"""

def _methods_inputs(wildcards):
    """Stats inputs for the methods report.

    Summarises the all-SNP (filtered), biallelic, and unlinked biallelic (thinned)
    datasets — counts and per-individual missingness for each.
    """
    proj = wildcards.project
    base = f"results/{proj}"
    return {
        "mainparams": "data/mainparams",
        "samples_to_keep": f"{base}/filtered_data/{proj}.samples_to_keep.txt",
        # Analysis (unlinked biallelic / thinned) dataset
        "vcf_stats": f"{base}/stats_vcf/thinned/{proj}.biallelic_snps.vcf_stats.txt",
        "depth_summary": f"{base}/stats_vcf/thinned/{proj}.biallelic_snps.depth_summary.txt",
        "imiss_thinned": f"{base}/stats_vcf/thinned/{proj}.biallelic_snps.imiss",
        # All-SNP (filtered) dataset
        "vcf_stats_filtered": f"{base}/stats_vcf/filtered/{proj}.filtered.vcf_stats.txt",
        "imiss_filtered": f"{base}/stats_vcf/filtered/{proj}.filtered.imiss",
        # Biallelic SNP dataset
        "vcf_stats_biallelic": f"{base}/stats_vcf/biallelic/{proj}.biallelic_snps.vcf_stats.txt",
        "imiss_biallelic": f"{base}/stats_vcf/biallelic/{proj}.biallelic_snps.imiss",
    }


rule generate_methods:
    input:
        unpack(_methods_inputs)
    output:
        methods = "results/{project}/methods.md"
    params:
        project    = lambda wildcards: wildcards.project,
        analyses   = lambda wildcards: config["projects"][wildcards.project].get("analyses", {}),
        parameters = lambda wildcards: config["projects"][wildcards.project].get("parameters", {}),
        envs_dir   = "workflow/envs",
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 5,
    script:
        "../scripts/generate_methods.py"
