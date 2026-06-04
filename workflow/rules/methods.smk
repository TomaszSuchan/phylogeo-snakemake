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

rule generate_methods:
    input:
        mainparams      = "data/mainparams",
        samples_to_keep = "results/{project}/filtered_data/{project}.samples_to_keep.txt",
        vcf_stats       = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.vcf_stats.txt",
        depth_summary   = "results/{project}/stats_vcf/thinned/{project}.biallelic_snps.depth_summary.txt",
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
