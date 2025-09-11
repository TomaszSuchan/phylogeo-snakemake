# Rule for PCAone analysis using EMU
rule pcaone_emu:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = "results/pcaone_EMU/PCA_EMU.eigvecs",
        pcaone_eigenvalues = "results/pcaone_EMU/PCA_EMU.eigvals"
    params:
        SVD_method = config.get("SVD_method", 3),
        output_prefix = "results/pcaone_EMU/PCA_EMU",
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
        mkdir -p $(dirname {params.output_prefix}) && \
        pcaone --threads {threads} \
        -d {params.SVD_method} \
        --emu \
        --bfile {params.bfile_prefix} \
        --out {params.output_prefix}
        """

# Rule to run PCAone analysis
rule pcaone:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = "results/pcaone/PCA.eigvecs",
        pcaone_eigenvalues = "results/pcaone/PCA.eigvals"
    params:
        SVD_method = config["PCA"].get("SVD_method", 3),
        output_prefix = "results/pcaone/PCA",
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
        mkdir -p $(dirname {params.output_prefix}) && \
        pcaone --threads {threads} \
        -d {params.SVD_method} \
        --bfile {params.bfile_prefix} \
        --out {params.output_prefix}
        """

# Rule to run PCAone for each mincov data threshold
rule pcaone_mincov:
    input:
        bed = "filtered_data/biallelic_snps_thinned_mincov{mincov}.bed",
        bim = "filtered_data/biallelic_snps_thinned_mincov{mincov}.bim",
        fam = "filtered_data/biallelic_snps_thinned_mincov{mincov}.fam"
    output:
        pcaone_eigenvectors = "results/pcaone_mincov{mincov}/PCA.eigvecs",
        pcaone_eigenvalues = "results/pcaone_mincov{mincov}/PCA.eigvals"
    params:
        SVD_method = config["PCA"].get("SVD_method", 3),
        output_prefix = lambda wildcards: f"results/pcaone_mincov{wildcards.mincov}/PCA",
        bfile_prefix = lambda wildcards: f"filtered_data/biallelic_snps_thinned_mincov{wildcards.mincov}"
    conda:
        "../envs/pcaone.yaml"
    threads: config["resources"]["pcaone"]["threads"]
    resources:
        mem_mb = config["resources"]["pcaone"]["mem_mb"],
        time = config["resources"]["pcaone"]["runtime"]
    shell:
        """
        mkdir -p $(dirname {params.output_prefix}) && \
        pcaone --threads {threads} \
        -d {params.SVD_method} \
        --bfile {params.bfile_prefix} \
        --out {params.output_prefix}
        """