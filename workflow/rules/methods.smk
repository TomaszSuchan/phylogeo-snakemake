"""
Rule to generate a manuscript-ready methods.md for each project.

The script reads the project config (analyses enabled/disabled, all parameters)
and the STRUCTURE mainparams file, then writes a single Markdown file with one
section per enabled analysis tool.  All config-derived parameter values are
filled in automatically; only runtime statistics (SNP counts, best K, software
versions) are left as [PLACEHOLDER].

Output: results/{project}/methods.md
"""

rule generate_methods:
    input:
        mainparams       = "data/mainparams",
        samples_to_keep  = "results/{project}/filtered_data/{project}.samples_to_keep.txt",
    output:
        methods = "results/{project}/methods.md"
    params:
        project    = lambda wildcards: wildcards.project,
        analyses   = lambda wildcards: config["projects"][wildcards.project].get("analyses", {}),
        parameters = lambda wildcards: config["projects"][wildcards.project].get("parameters", {}),
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 5,
    script:
        "../scripts/generate_methods.py"
