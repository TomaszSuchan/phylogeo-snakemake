# Rule for PCAone analysis using EMU
rule pcaone_emu:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = "results/{project}/pcaone_EMU/{project}.PCA_EMU.eigvecs",
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
        eigenvalues = "results/{project}/pcaone_miss{miss}/{project}.PCA_miss{miss}.eigvals"
    log:
        "logs/{project}/pcaone_miss_{miss}.log"
    benchmark:
        "benchmarks/{project}/pcaone_miss_{miss}.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("SVD_method", 3),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = lambda wildcards: f"results/{wildcards.project}/pcaone_miss{wildcards.miss}/{wildcards.project}.PCA",
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

# Rule to plot PCA
rule plot_pca:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        "results/{project}/pcaone/plots/{project}.PCA-{color_by}.pdf"
    params:
        pc1 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pca_plot"]["pc1"],
        pc2 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pca_plot"]["pc2"],
        color_by = lambda wildcards: wildcards.color_by
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca.R"