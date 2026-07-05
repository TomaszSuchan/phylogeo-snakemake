# Align ancestry clusters (Q-matrix columns) so a given biological cluster keeps the
# same column - and therefore the same colour - across K within a method and across
# methods. See workflow/scripts/align_clusters.R for the algorithm.
#
# A single generic rule serves all structure-like methods via the {method} path
# wildcard. Alignment needs every K at once, so the rule takes all K as input and
# writes the aligned matrix for the requested K; Snakemake runs it once per K. Methods
# other than the reference additionally take the reference method's aligned matrices
# as the anchor (wired by _align_inputs / _ancestry_reference_method in the Snakefile).

rule align_ancestry_clusters:
    input:
        unpack(_align_inputs)
    output:
        aligned = "results/{project}/{method}/aligned/{project}.{method}.K{k}.Qmatrix.aligned.txt"
    wildcard_constraints:
        method = "|".join(ANCESTRY_METHODS),
        k = PLOT_K_WILDCARD_CONSTRAINT
    params:
        method = lambda wildcards: wildcards.method,
        k_values = lambda wildcards: _ancestry_align_kvals(wildcards.project),
        ref_k_values = lambda wildcards: (
            _ancestry_align_kvals(wildcards.project)
            if _ancestry_reference_method(wildcards.project) != wildcards.method
            else []
        ),
        target_k = lambda wildcards: wildcards.k
    log:
        "logs/{project}/align_clusters.{method}.K{k}.log"
    conda:
        "../envs/r-cluster-align.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    script:
        "../scripts/align_clusters.R"
