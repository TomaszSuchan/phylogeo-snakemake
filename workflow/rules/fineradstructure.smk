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
        "logs/{project}/prepare_finestr_input.log"
    benchmark:
        "benchmarks/{project}/prepare_finestr_input.txt"
    conda:
        "../envs/fineradstructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        RADpainter hapsFromVCF {input.vcf} > {output.finestr_input} &> {log}
        """

rule fineradstructure_paint:
    input:
        finestr_input = rules.fineradstructure_prepare_input.output.finestr_input
    output:
        finestr_chunks = "results/{project}/fineradstructure/{project}_chunks.out"
    log:
        "logs/{project}/fineradstructure/run_fineradstructure.log"
    benchmark:
        "benchmarks/{project}/run_fineradstructure.txt"
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
        "logs/{project}/fineradstructure/cluster_fineradstructure.log"
    benchmark:
        "benchmarks/{project}/cluster_fineradstructure.txt"
    conda:
        "../envs/fineradstructure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        finestructure -x {params.mcmc_iterations} -y {params.burnin} -z {params.thinning} -X -Y {input.finestr_chunks} {output.finestr_mcmc} &> {log}
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
        "logs/{project}/fineradstructure/fineradstructure_tree.log"
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
        finestructure -m T -x {params.mcmc_iterations} {input.finestr_chunks} {input.finestr_mcmc} {output.finestr_mcmcTree} &> {log}
        """

rule fineradstructure_plot:
    """
    Plot fineRADstructure results - visualizes the co-ancestry tree from MCMC output
    """
    input:
        mcmcTree = rules.fineradstructure_tree.output.finestr_mcmcTree,
        chunks = rules.fineradstructure_paint.output.finestr_chunks
    output:
        pdf = "results/{project}/fineradstructure/plots/{project}.fineradstructure_tree.pdf",
        rds = "results/{project}/fineradstructure/plots/{project}.fineradstructure_tree.rds"
    log:
        "logs/{project}/fineradstructure/fineradstructure_plot.log"
    benchmark:
        "benchmarks/{project}/fineradstructure_plot.txt"
    conda:
        "../envs/r-ggtree.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/fineradstructure_plot.R"