# Rule for ADMIXTURE
rule admixture:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        Q = "results/{project}/admixture/{project}.biallelic_snps_thinned.{k}.Q",
        P = "results/{project}/admixture/{project}.biallelic_snps_thinned.{k}.P"
    log:
        "logs/{project}/admixture.{k}.log"
    benchmark:
        "benchmarks/{project}/admixture.K{k}.txt"
    params:
        output_dir = "results/{project}/admixture/",
        prefix = lambda wildcards: f"{wildcards.project}.biallelic_snps_thinned"
    conda:
        "../envs/admixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["admixture"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["admixture"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["admixture"]["runtime"]
    shell:
        """
        mkdir -p {params.output_dir}
        admixture --cv -j{threads} {input.bed} {wildcards.k} > {log} 2>&1
        mv {params.prefix}.{wildcards.k}.Q {params.output_dir}/
        mv {params.prefix}.{wildcards.k}.P {params.output_dir}/
        """

# Rule to choose optimal K for admixture
rule admixture_chooseK:
    input:
        lambda wildcards: expand(
            "logs/{project}/admixture.{k}.log",
            project=wildcards.project,
            k=config["projects"][wildcards.project]["parameters"]["k_values"]
        )
    output:
        choose_k_results = "results/{project}/admixture/{project}.admixture.chooseK_results.txt",
        cv_summary = "results/{project}/admixture/{project}.admixture.cv_summary.txt"
    log:
        "logs/{project}/admixture_choose_k.log"
    benchmark:
        "benchmarks/{project}/admixture_choose_k.txt"
    script:
        "../scripts/admixture_choose_k.py"


rule plot_admixture_cv:
    """
    Plot ADMIXTURE cross-validation error across K (evanno-style ggplot aesthetics).
    """
    input:
        cv_summary = rules.admixture_chooseK.output.cv_summary
    output:
        pdf = "results/{project}/admixture/plots/{project}.admixture.cv_plot.pdf",
        rds = "results/{project}/admixture/plots/{project}.admixture.cv_plot.rds"
    log:
        "logs/{project}/plot_admixture_cv.log"
    params:
        ylab = "Cross-validation error",
        width = lambda wildcards: _choose_k_plot_param(wildcards, "width", 25.4),
        height = lambda wildcards: _choose_k_plot_param(wildcards, "height", 12.7),
        dpi = lambda wildcards: _choose_k_plot_param(wildcards, "dpi", 300)
    conda:
        "../envs/mapmixture.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_choose_k_score.R"
