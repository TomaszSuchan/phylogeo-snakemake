# DAPC (Discriminant Analysis of Principal Components) rule - run for each K
rule dapc_analysis:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    wildcard_constraints:
        # K>=2 as decimal strings: [2-9]\d* wrongly rejects 10,11,... (first digit is 1).
        k = r"(?:[2-9]|[1-9]\d+)"
    output:
        results_rds = "results/{project}/dapc/{project}.dapc.K{k}.results.rds",
        cluster_assignments = "results/{project}/dapc/{project}.dapc.K{k}.cluster_assignments.txt",
        membership_probs = "results/{project}/dapc/{project}.dapc.K{k}.membership_probs.txt",
        scatter_plot = "results/{project}/dapc/plots/{project}.dapc.K{k}.scatter.pdf",
        scatter_plot_rds = "results/{project}/dapc/plots/{project}.dapc.K{k}.scatter.rds"
    params:
        k = lambda wildcards: int(wildcards.k),
        n_pca = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("dapc", {}).get("n_pca", "retained"),
        n_da = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("dapc", {}).get("n_da", "all"),
        criterion = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("dapc", {}).get("criterion", "diffNgroup")
    conda:
        "../envs/adegenet.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["dapc"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["dapc"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["dapc"]["runtime"]
    log:
        "logs/{project}/dapc.K{k}.log.txt"
    benchmark:
        "benchmarks/{project}/dapc.K{k}.txt"
    script:
        "../scripts/dapc_analysis.R"

# Rule to plot BIC values parsed from per-K DAPC logs
rule dapc_bic_plot_from_log:
    input:
        log_files = lambda wildcards: expand(
            "logs/{project}/dapc.K{k}.log.txt",
            project=wildcards.project,
            k=[k for k in config["projects"][wildcards.project]["parameters"]["k_values"] if int(k) >= 2]
        )
    output:
        plot = "results/{project}/dapc/plots/{project}.dapc.criterion_plot.pdf",
        plot_rds = "results/{project}/dapc/plots/{project}.dapc.criterion_plot.rds"
    conda:
        "../envs/r-plot.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    log:
        "logs/{project}/dapc_bic_plot_from_log.log"
    params:
        width = lambda wildcards: _choose_k_plot_param(wildcards, "width", 10),
        height = lambda wildcards: _choose_k_plot_param(wildcards, "height", 5),
        dpi = lambda wildcards: _choose_k_plot_param(wildcards, "dpi", 300)
    script:
        "../scripts/plot_dapc_bic_from_log.R"

# Rule to plot DAPC results on map with mapmixture (for specific K)
rule mapmixture_dapc:
    input:
        unpack(_dapc_map_inputs),
    output:
        plot = "results/{project}/dapc/plots/{project}.dapc.K{k}.map.pdf",
        plot_rds = "results/{project}/dapc/plots/{project}.dapc.K{k}.map.rds"
    params:
        unpack(lambda wildcards: _mapmixture_map_params(wildcards))
    log:
        "logs/{project}/mapmixture_dapc.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_dapc_map.R"
