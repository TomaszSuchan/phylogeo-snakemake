# Rule for PCAone analysis using EMU
rule pcaone_emu:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = "pcaone_EMU/PCA_EMU.eigvecs",
        pcaone_eigenvalues = "pcaone_EMU/PCA_EMU.eigvals"
    params:
        SVD_method = config.get("SVD_method", 3),
        output_prefix = "pcaone_EMU/PCA_EMU",
        # Get the bfile prefix (remove .bed extension)
        bfile_prefix = lambda wildcards, input: input.bed.replace('.bed', '')
    conda:
        "../envs/pcaone.yaml"
    threads: 4
    resources:
        mem_mb = 8000,
        time = "1:00:00"
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
        pcaone_eigenvectors = "pcaone/PCA.eigvecs",
        pcaone_eigenvalues = "pcaone/PCA.eigvals"
    params:
        SVD_method = config.get("SVD_method", 3),
        output_prefix = "pcaone/PCA",
        # Get the bfile prefix (remove .bed extension)
        bfile_prefix = lambda wildcards, input: input.bed.replace('.bed', '')
    conda:
        "../envs/pcaone.yaml"
    threads: 4
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

# Rule to run PCAone analyses
rule run_PCAone:
    input:
        rules.pcaone_emu.output,
        rules.pcaone.output