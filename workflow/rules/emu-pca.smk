# EMU-PCA: EM-PCA for ultra-low coverage sequencing data (Rosemeis/emu)
# Uses PLINK binary files; outputs .eigvecs (tab-sep with #FID IID PC1...) and .eigvals
# https://github.com/Rosemeis/emu

rule emu_pca:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        emupca_eigenvectors = "results/{project}/emupca/{project}.emupca.eigvecs",
        emupca_eigenvalues = "results/{project}/emupca/{project}.emupca.eigvals"
    log:
        "logs/{project}/emu_pca.log"
    benchmark:
        "benchmarks/{project}/emu_pca.txt"
    params:
        eig = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("EMU", {}).get("eig", 2),
        eig_out = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("EMU", {}).get("eig_out", 10),
        output_prefix = "results/{project}/emupca/{project}.emupca",
        bfile_prefix = lambda wildcards, input: input.bed.replace(".bed", "")
    conda:
        "../envs/emu-pca.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["emupca"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["emupca"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["emupca"]["runtime"]
    shell:
        """
        emu --bfile {params.bfile_prefix} \
            --eig {params.eig} \
            --eig-out {params.eig_out} \
            --threads {threads} \
            --out {params.output_prefix} &> {log}
        """

# Rule to plot EMU-PCA (colored by population/metadata)
rule plot_emupca_colored:
    input:
        eigvecs = rules.emu_pca.output.emupca_eigenvectors,
        eigvals = rules.emu_pca.output.emupca_eigenvalues,
        indpopdata = rules.generate_popdata.output.indpopdata,
        indmiss = rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf = "results/{project}/emupca/plots/{project}.emupca-PC{pc1}_PC{pc2}-{color_by}.pdf",
        rds = "results/{project}/emupca/plots/{project}.emupca-PC{pc1}_PC{pc2}-{color_by}.rds"
    log:
        "logs/{project}/plot_emupca_colored_PC{pc1}_PC{pc2}_{color_by}.log"
    wildcard_constraints:
        color_by = "(?!labeled|missing).*"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        color_by = lambda wildcards: wildcards.color_by,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        plot_type = "colored"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_emupca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot EMU-PCA with labels only
rule plot_emupca_labeled:
    input:
        eigvecs = rules.emu_pca.output.emupca_eigenvectors,
        eigvals = rules.emu_pca.output.emupca_eigenvalues,
        indpopdata = rules.generate_popdata.output.indpopdata,
        indmiss = rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf = "results/{project}/emupca/plots/{project}.emupca-PC{pc1}_PC{pc2}-labeled.pdf",
        rds = "results/{project}/emupca/plots/{project}.emupca-PC{pc1}_PC{pc2}-labeled.rds"
    log:
        "logs/{project}/plot_emupca_labeled_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "labeled"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_emupca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot EMU-PCA colored by missing data
rule plot_emupca_missing:
    input:
        eigvecs = rules.emu_pca.output.emupca_eigenvectors,
        eigvals = rules.emu_pca.output.emupca_eigenvalues,
        indpopdata = rules.generate_popdata.output.indpopdata,
        indmiss = rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf = "results/{project}/emupca/plots/{project}.emupca-PC{pc1}_PC{pc2}-missing.pdf",
        rds = "results/{project}/emupca/plots/{project}.emupca-PC{pc1}_PC{pc2}-missing.rds"
    log:
        "logs/{project}/plot_emupca_missing_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "missing"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_emupca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot EMU-PCA facet (all PC combinations, colored)
rule plot_emupca_facet_colored:
    input:
        eigvecs = rules.emu_pca.output.emupca_eigenvectors,
        eigvals = rules.emu_pca.output.emupca_eigenvalues,
        indpopdata = rules.generate_popdata.output.indpopdata,
        indmiss = rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf = "results/{project}/emupca/plots/{project}.emupca-facet-{color_by}.pdf",
        rds = "results/{project}/emupca/plots/{project}.emupca-facet-{color_by}.rds"
    log:
        "logs/{project}/plot_emupca_facet_colored_{color_by}.log"
    wildcard_constraints:
        color_by = "(?!labeled|missing).*"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "colored"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_emupca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot EMU-PCA facet (all PC combinations, labeled)
rule plot_emupca_facet_labeled:
    input:
        eigvecs = rules.emu_pca.output.emupca_eigenvectors,
        eigvals = rules.emu_pca.output.emupca_eigenvalues,
        indpopdata = rules.generate_popdata.output.indpopdata,
        indmiss = rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf = "results/{project}/emupca/plots/{project}.emupca-facet-labeled.pdf",
        rds = "results/{project}/emupca/plots/{project}.emupca-facet-labeled.rds"
    log:
        "logs/{project}/plot_emupca_facet_labeled.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "labeled"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_emupca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot EMU-PCA facet (all PC combinations, missing)
rule plot_emupca_facet_missing:
    input:
        eigvecs = rules.emu_pca.output.emupca_eigenvectors,
        eigvals = rules.emu_pca.output.emupca_eigenvalues,
        indpopdata = rules.generate_popdata.output.indpopdata,
        indmiss = rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf = "results/{project}/emupca/plots/{project}.emupca-facet-missing.pdf",
        rds = "results/{project}/emupca/plots/{project}.emupca-facet-missing.rds"
    log:
        "logs/{project}/plot_emupca_facet_missing.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "missing"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_emupca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"
