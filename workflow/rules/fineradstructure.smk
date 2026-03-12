## fineRADstructure analysis targets
rule fineradstructure_all:
    """Target rule to run complete fineRADstructure analysis"""
    input:
        finestr_mcmcTree = "results/{project}/fineradstructure/{project}_chunks.mcmcTree.xml"
    default_target: True

rule fineradstructure_prepare_input:
    input:
        vcf = rules.select_biallelic_snps.output.biallelic_vcf
    output:
        finestr_input = "results/{project}/fineradstructure/{project}.input"
    log:
        "logs/{project}/fineradstructure_prepare_input.log"
    benchmark:
        "benchmarks/{project}/fineradstructure_prepare_input.txt"
    conda:
        "../envs/fineradstructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        RADpainter hapsFromVCF {input.vcf} > {output.finestr_input} 2> {log}
        """

rule fineradstructure_paint:
    input:
        finestr_input = rules.fineradstructure_prepare_input.output.finestr_input
    output:
        finestr_chunks = "results/{project}/fineradstructure/{project}_chunks.out"
    log:
        "logs/{project}/fineradstructure_paint.log"
    benchmark:
        "benchmarks/{project}/fineradstructure_paint.txt"
    conda:
        "../envs/fineradstructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        RADpainter paint {input.finestr_input} &> {log}
        """

rule fineradstructure_cluster:
    input:
        finestr_chunks = rules.fineradstructure_paint.output.finestr_chunks
    output:
        finestr_mcmc = "results/{project}/fineradstructure/{project}_chunks.mcmc.xml"
    params:
        mcmc_iterations = lambda wildcards: config["projects"][wildcards.project]["parameters"]["fineradstructure"]["cluster"]["mcmc_iterations"],
        burnin = lambda wildcards: config["projects"][wildcards.project]["parameters"]["fineradstructure"]["cluster"]["burnin"],
        thinning = lambda wildcards: config["projects"][wildcards.project]["parameters"]["fineradstructure"]["cluster"]["thinning"]
    log:
        "logs/{project}/fineradstructure_cluster.log"
    benchmark:
        "benchmarks/{project}/fineradstructure_cluster.txt"
    conda:
        "../envs/fineradstructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["fineradstructure_cluster"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["fineradstructure_cluster"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["fineradstructure_cluster"]["runtime"]
    shell:
        """
        finestructure -Y -x {params.mcmc_iterations} -y {params.burnin} -z {params.thinning} -X -Y {input.finestr_chunks} {output.finestr_mcmc} &> {log}
        """

rule fineradstructure_tree:
    input:
        finestr_chunks = rules.fineradstructure_paint.output.finestr_chunks,
        finestr_mcmc = rules.fineradstructure_cluster.output.finestr_mcmc
    output:
        finestr_mcmcTree = "results/{project}/fineradstructure/{project}_chunks.mcmcTree.xml"
    params:
        mcmc_iterations = lambda wildcards: config["projects"][wildcards.project]["parameters"]["fineradstructure"]["tree"]["mcmc_iterations"]
    log:
        "logs/{project}/fineradstructure_tree.log"
    benchmark:
        "benchmarks/{project}/fineradstructure_tree.txt"
    conda:
        "../envs/fineradstructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        finestructure -Y -m T -x {params.mcmc_iterations} {input.finestr_chunks} {input.finestr_mcmc} {output.finestr_mcmcTree} &> {log}
        """

rule fineradstructure_plot:
    """
    Plot fineRADstructure results using the upstream-style heatmap workflow
    """
    input:
        mcmcTree = rules.fineradstructure_tree.output.finestr_mcmcTree,
        mcmc = rules.fineradstructure_cluster.output.finestr_mcmc,
        chunks = rules.fineradstructure_paint.output.finestr_chunks,
        indpopdata = rules.generate_popdata.output.indpopdata
    output:
        simple_pdf = "results/{project}/fineradstructure/plots/{project}.SimpleCoancestry.pdf",
        popavg_pdf = "results/{project}/fineradstructure/plots/{project}.PopAveragedCoancestry.pdf",
        labeled_pdf = "results/{project}/fineradstructure/plots/{project}.PopAveragedCoancestryLabeled.pdf",
        rds = "results/{project}/fineradstructure/plots/{project}.fineradstructure_plots.rds"
    params:
        label_by = lambda wildcards: config["projects"][wildcards.project]["parameters"]["fineradstructure"].get("plot", {}).get("label_by", "Site"),
        max_indv = lambda wildcards: config["projects"][wildcards.project]["parameters"]["fineradstructure"].get("plot", {}).get("max_indv", 10000),
        max_pop = lambda wildcards: config["projects"][wildcards.project]["parameters"]["fineradstructure"].get("plot", {}).get("max_pop", 10000),
        max_label_values = lambda wildcards: config["projects"][wildcards.project]["parameters"]["fineradstructure"].get("plot", {}).get("max_label_values", 3)
    log:
        "logs/{project}/fineradstructure_plot.log"
    benchmark:
        "benchmarks/{project}/fineradstructure_plot.txt"
    conda:
        "../envs/fineradstructure_plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/fineradstructure_plot.R"