# Rule for PCAone analysis using EMU
rule pcaone_emu:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = "results/{project}/pcaone_EMU/{project}.PCA_EMU.eigvecs",
        pcaone_eigenvectors2 = "results/{project}/pcaone_EMU/{project}.PCA_EMU.eigvecs2",
        pcaone_eigenvalues = "results/{project}/pcaone_EMU/{project}.PCA_EMU.eigvals"
    log:
        "logs/{project}/pcaone_emu.log"
    benchmark:
        "benchmarks/{project}/pcaone_emu.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("SVD_method", 3),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = "results/{project}/pcaone_EMU/{project}.PCA_EMU",
        # Get the bfile prefix (remove .bed extension)
        bfile_prefix = lambda wildcards, input: input.bed.replace('.bed', '')
    conda:
        "../envs/pcaone.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone_emu"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone_emu"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone_emu"]["runtime"]
    shell:
        """
        PCAone --threads {threads} \
        -d {params.SVD_method} \
        --pc {params.PCnum} \
        --emu \
        --bfile {params.bfile_prefix} \
        --out {params.output_prefix} &> {log}
        """

# Rule to run PCAone analysis
rule pcaone:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = "results/{project}/pcaone/{project}.PCA.eigvecs",
        pcaone_eigenvectors2 = "results/{project}/pcaone/{project}.PCA.eigvecs2",
        pcaone_eigenvalues = "results/{project}/pcaone/{project}.PCA.eigvals"
    log:
        "logs/{project}/pcaone.log"
    benchmark:
        "benchmarks/{project}/pcaone.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("SVD_method", 3),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = "results/{project}/pcaone/{project}.PCA",
        # Get the bfile prefix (remove .bed extension)
        bfile_prefix = lambda wildcards, input: input.bed.replace('.bed', '')
    conda:
        "../envs/pcaone.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["runtime"]
    shell:
        """
        PCAone --threads {threads} \
        -d {params.SVD_method} \
        --pc {params.PCnum} \
        --bfile {params.bfile_prefix} \
        --out {params.output_prefix} &> {log}
        """

# Rule to run PCAone for each miss data threshold
rule pcaone_miss:
    input:
        bed = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.bed",
        bim = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.bim",
        fam = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.fam"
    output:
        eigenvectors = "results/{project}/pcaone_miss{miss}/{project}.PCA_miss{miss}.eigvecs",
        eigenvectors2 = "results/{project}/pcaone_miss{miss}/{project}.PCA_miss{miss}.eigvecs2",
        eigenvalues = "results/{project}/pcaone_miss{miss}/{project}.PCA_miss{miss}.eigvals"
    log:
        "logs/{project}/pcaone_miss_{miss}.log"
    benchmark:
        "benchmarks/{project}/pcaone_miss_{miss}.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("SVD_method", 3),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = lambda wildcards: f"results/{wildcards.project}/pcaone_miss{wildcards.miss}/{wildcards.project}.PCA_miss{wildcards.miss}",
        bfile_prefix = lambda wildcards: f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_thinned_miss{wildcards.miss}"
    conda:
        "../envs/pcaone.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["runtime"]
    shell:
        """
        PCAone --threads {threads} \
        -d {params.SVD_method} \
        --pc {params.PCnum} \
        --bfile {params.bfile_prefix} \
        --out {params.output_prefix} &> {log}
        """

# Rule to plot PCA (colored by population/metadata)
rule plot_pca_colored:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-{color_by}.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_colored_PC{pc1}_PC{pc2}_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
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
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA with labels only
rule plot_pca_labeled:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-labeled.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-labeled.rds"
    log:
        "logs/{project}/plot_pca_labeled_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "labeled"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA colored by missing data
rule plot_pca_missing:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-missing.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-missing.rds"
    log:
        "logs/{project}/plot_pca_missing_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "missing"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA facet (all PC combinations, colored)
rule plot_pca_facet_colored:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-facet-{color_by}.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-facet-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_facet_colored_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        pca_colors = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pca_colors", None),
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "colored"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA facet (all PC combinations, labeled)
rule plot_pca_facet_labeled:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-facet-labeled.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-facet-labeled.rds"
    log:
        "logs/{project}/plot_pca_facet_labeled.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "labeled"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA facet (all PC combinations, missing)
rule plot_pca_facet_missing:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-facet-missing.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-facet-missing.rds"
    log:
        "logs/{project}/plot_pca_facet_missing.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "missing"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA EMU (colored by population/metadata)
rule plot_pca_emu_colored:
    input:
        eigvecs=rules.pcaone_emu.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone_emu.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    output:
        pdf="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-{color_by}.pdf",
        rds="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_emu_colored_PC{pc1}_PC{pc2}_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
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
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA EMU with labels only
rule plot_pca_emu_labeled:
    input:
        eigvecs=rules.pcaone_emu.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone_emu.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    output:
        pdf="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-labeled.pdf",
        rds="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-labeled.rds"
    log:
        "logs/{project}/plot_pca_emu_labeled_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "labeled"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA EMU colored by missing data
rule plot_pca_emu_missing:
    input:
        eigvecs=rules.pcaone_emu.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone_emu.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    output:
        pdf="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-missing.pdf",
        rds="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-missing.rds"
    log:
        "logs/{project}/plot_pca_emu_missing_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "missing"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA for each miss threshold (colored by population/metadata)
rule plot_pca_miss_colored:
    input:
        eigvecs=rules.pcaone_miss.output.eigenvectors2,
        eigvals=rules.pcaone_miss.output.eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.imiss"
    output:
        pdf="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-{color_by}.pdf",
        rds="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_miss{miss}_colored_PC{pc1}_PC{pc2}_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
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
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA miss with labels only
rule plot_pca_miss_labeled:
    input:
        eigvecs=rules.pcaone_miss.output.eigenvectors2,
        eigvals=rules.pcaone_miss.output.eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.imiss"
    output:
        pdf="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-labeled.pdf",
        rds="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-labeled.rds"
    log:
        "logs/{project}/plot_pca_miss{miss}_labeled_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "labeled"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA miss colored by missing data
rule plot_pca_miss_missing:
    input:
        eigvecs=rules.pcaone_miss.output.eigenvectors2,
        eigvals=rules.pcaone_miss.output.eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss="results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.imiss"
    output:
        pdf="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-missing.pdf",
        rds="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-missing.rds"
    log:
        "logs/{project}/plot_pca_miss{miss}_missing_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "missing"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"