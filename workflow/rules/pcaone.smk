# Rule for PCAone analysis using EMU
rule pcaone_emu:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = "{project}/pcaone_EMU/PCA_EMU.eigvecs",
        pcaone_eigenvalues = "{project}/pcaone_EMU/PCA_EMU.eigvals"
    log:
        "{project}/logs/pcaone_emu.log"
    benchmark:
        "{project}/benchmarks/pcaone_emu.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("SVD_method", 3),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = "{project}/pcaone_EMU/PCA_EMU",
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
        pcaone_eigenvectors = "{project}/pcaone/PCA.eigvecs",
        pcaone_eigenvectors2 = "{project}/pcaone/PCA.eigvecs2",
        pcaone_eigenvalues = "{project}/pcaone/PCA.eigvals"
    log:
        "{project}/logs/pcaone.log"
    benchmark:
        "{project}/benchmarks/pcaone.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("SVD_method", 3),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = "{project}/pcaone/PCA",
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
        bed = "{project}/filtered_data/biallelic_snps_thinned_miss{miss}.bed",
        bim = "{project}/filtered_data/biallelic_snps_thinned_miss{miss}.bim",
        fam = "{project}/filtered_data/biallelic_snps_thinned_miss{miss}.fam"
    output:
        eigenvectors = "{project}/pcaone_miss{miss}/PCA_miss{miss}.eigvecs",
        eigenvalues = "{project}/pcaone_miss{miss}/PCA_miss{miss}.eigvals"
    log:
        "{project}/logs/pcaone_miss_{miss}.log"
    benchmark:
        "{project}/benchmarks/pcaone_miss_{miss}.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("SVD_method", 3),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = lambda wildcards: f"{config['analysis_name']}/pcaone_miss{wildcards.miss}/PCA",
        bfile_prefix = lambda wildcards: f"{config['analysis_name']}/filtered_data/biallelic_snps_thinned_miss{wildcards.miss}"
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

# Rule to plot PCA
rule plot_pca:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        "{project}/final_plots/PCA-{color_by}.pdf"
    params:
        pc1 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pca_plot"]["pc1"],
        pc2 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pca_plot"]["pc2"],
        color_by = lambda wildcards: wildcards.color_by
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca.R"