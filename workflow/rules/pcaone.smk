# Rule for PCAone analysis using EMU
rule pcaone_emu:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = config["analysis_name"] + "/pcaone_EMU/PCA_EMU.eigvecs",
        pcaone_eigenvalues = config["analysis_name"] + "/pcaone_EMU/PCA_EMU.eigvals"
    log:
        config["analysis_name"] + "/logs/pcaone_emu.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/pcaone_emu.txt"
    params:
        SVD_method = config["PCAone"].get("SVD_method", 3),
        PCnum = config["PCAone"].get("PCnum", 10),
        output_prefix = config["analysis_name"] + "/pcaone_EMU/PCA_EMU",
        # Get the bfile prefix (remove .bed extension)
        bfile_prefix = lambda wildcards, input: input.bed.replace('.bed', '')
    conda:
        "../envs/pcaone.yaml"
    threads: config["resources"]["pcaone_emu"]["threads"]
    resources:
        mem_mb = config["resources"]["pcaone_emu"]["mem_mb"],
        time = config["resources"]["pcaone_emu"]["runtime"]
    threads: 12
    resources:
        mem_mb = 8000,
        time = "72:00:00"
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
        pcaone_eigenvectors = config["analysis_name"] + "/pcaone/PCA.eigvecs",
        pcaone_eigenvectors2 = config["analysis_name"] + "/pcaone/PCA.eigvecs2", 
        pcaone_eigenvalues = config["analysis_name"] + "/pcaone/PCA.eigvals"
    log:
        config["analysis_name"] + "/logs/pcaone.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/pcaone.txt"
    params:
        SVD_method = config["PCAone"].get("SVD_method", 3),
        PCnum = config["PCAone"].get("PCnum", 10),
        output_prefix = config["analysis_name"] + "/pcaone/PCA",
        # Get the bfile prefix (remove .bed extension)
        bfile_prefix = lambda wildcards, input: input.bed.replace('.bed', '')
    conda:
        "../envs/pcaone.yaml"
    threads: config["resources"]["pcaone"]["threads"]
    resources:
        mem_mb = config["resources"]["pcaone"]["mem_mb"],
        time = config["resources"]["pcaone"]["runtime"]
    threads: 1
    resources:
        mem_mb = 8000,
        time = "1:00:00"
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
        bed = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned_miss{miss}.bed",
        bim = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned_miss{miss}.bim",
        fam = config["analysis_name"] + "/filtered_data/biallelic_snps_thinned_miss{miss}.fam"
    output:
        eigenvectors = config["analysis_name"] + "/pcaone_miss{miss}/PCA_miss{miss}.eigvecs",
        eigenvalues = config["analysis_name"] + "/pcaone_miss{miss}/PCA_miss{miss}.eigvals"
    log:
        config["analysis_name"] + "/logs/pcaone_miss_{miss}.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/pcaone_miss_{miss}.txt"
    params:
        SVD_method = config["PCAone"].get("SVD_method", 3),
        PCnum = config["PCAone"].get("PCnum", 10),
        output_prefix = lambda wildcards: f"{config['analysis_name']}/pcaone_miss{wildcards.miss}/PCA",
        bfile_prefix = lambda wildcards: f"{config['analysis_name']}/filtered_data/biallelic_snps_thinned_miss{wildcards.miss}"
    conda:
        "../envs/pcaone.yaml"
    threads: config["resources"]["pcaone"]["threads"]
    resources:
        mem_mb = config["resources"]["pcaone"]["mem_mb"],
        time = config["resources"]["pcaone"]["runtime"]
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
        popdata=rules.generate_popdata.output.popdata
    output:
        "{analysis}/final_plots/PCA-{color_by}.pdf"
    params:
        pc1=config["pca_plot"]["pc1"],
        pc2=config["pca_plot"]["pc2"],
        color_by=config["pca_plot"]["color_by"]
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca.R"
