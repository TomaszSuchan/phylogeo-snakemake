# Rules for plotting STRUCTURE, fastStructure, and ADMIXTURE results on maps using mapmixture

# Rule to install mapmixture once (all plotting rules depend on this)
rule install_mapmixture:
    output:
        touch(".snakemake/mapmixture_installed")
    conda:
        "../envs/mapmixture.yaml"
    shell:
        """
        Rscript -e 'lib <- .libPaths()[1]; options(repos="https://cloud.r-project.org/"); pkgs <- c("mapmixture", "elevatr", "rnaturalearth", "RColorBrewer"); for (p in pkgs) if (!requireNamespace(p, quietly=TRUE, lib.loc=lib)) install.packages(p, lib=lib)'
        """


# Rule to plot STRUCTURE results on map with mapmixture
rule mapmixture_structure:
    input:
        unpack(_structure_map_inputs),
    output:
        plot = "results/{project}/structure/plots/{project}.structure.K{k}.map.pdf",
        plot_rds = "results/{project}/structure/plots/{project}.structure.K{k}.map.rds"
    params:
        lambda wildcards: _mapmixture_map_rule_params(wildcards, "structure", "structure")
    log:
        "logs/{project}/mapmixture_structure.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"

# Rule to plot fastStructure results on map with mapmixture
rule mapmixture_faststructure:
    input:
        unpack(_faststructure_map_inputs),
    output:
        plot = "results/{project}/faststructure/plots/{project}.faststructure.K{k}.map.pdf",
        plot_rds = "results/{project}/faststructure/plots/{project}.faststructure.K{k}.map.rds"
    params:
        lambda wildcards: _mapmixture_map_rule_params(wildcards, "faststructure", "faststructure")
    log:
        "logs/{project}/mapmixture_faststructure.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"

# Rule to plot ADMIXTURE results on map with mapmixture
rule mapmixture_admixture:
    input:
        unpack(_admixture_map_inputs),
    output:
        plot = "results/{project}/admixture/plots/{project}.admixture.K{k}.map.pdf",
        plot_rds = "results/{project}/admixture/plots/{project}.admixture.K{k}.map.rds"
    params:
        lambda wildcards: _mapmixture_map_rule_params(wildcards, "admixture", "admixture")
    log:
        "logs/{project}/mapmixture_admixture.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_map.R"

# Rule to generate structure barplots from STRUCTURE results
rule barplot_structure:
    input:
        qmatrix = "results/{project}/structure/{project}.structure.K{k}.Qmatrix.txt",
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_mapmixture.output
    output:
        barplot = "results/{project}/structure/plots/{project}.structure.K{k}.barplot.pdf",
        barplot_rds = "results/{project}/structure/plots/{project}.structure.K{k}.barplot.rds"
    params:
        lambda wildcards: _barplot_rule_params(wildcards, "structure", "structure")
    log:
        "logs/{project}/structure_barplot.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_barplots.R"

# Rule to generate structure barplots from fastStructure results
rule barplot_faststructure:
    input:
        qmatrix = "results/{project}/faststructure/{project}.faststructure.{k}.meanQ",
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_mapmixture.output
    output:
        barplot = "results/{project}/faststructure/plots/{project}.faststructure.K{k}.barplot.pdf",
        barplot_rds = "results/{project}/faststructure/plots/{project}.faststructure.K{k}.barplot.rds"
    params:
        lambda wildcards: _barplot_rule_params(wildcards, "faststructure", "faststructure")
    log:
        "logs/{project}/faststructure_barplot.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_barplots.R"

# Rule to generate structure barplots from ADMIXTURE results
rule barplot_admixture:
    input:
        qmatrix = "results/{project}/admixture/{project}.biallelic_snps_thinned.{k}.Q",
        indpopdata = rules.generate_popdata.output.indpopdata,
        install = rules.install_mapmixture.output
    output:
        barplot = "results/{project}/admixture/plots/{project}.admixture.K{k}.barplot.pdf",
        barplot_rds = "results/{project}/admixture/plots/{project}.admixture.K{k}.barplot.rds"
    params:
        lambda wildcards: _barplot_rule_params(wildcards, "admixture", "admixture")
    log:
        "logs/{project}/admixture_barplot.K{k}.log"
    conda:
        "../envs/mapmixture.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_structure_barplots.R"
