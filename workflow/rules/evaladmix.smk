"""
Rules for evalAdmix - evaluation of admixture analysis results.
evalAdmix evaluates the fit of admixture models by calculating correlation of residuals.
Works with ADMIXTURE, fastStructure, and STRUCTURE outputs.
"""

# Rule to download visFuns.R from evalAdmix repository
rule download_visfuns:
    output:
        visfuns = "workflow/scripts/visFuns.R"
    shell:
        """
        # Create directory if it doesn't exist
        mkdir -p $(dirname {output.visfuns})
        
        # Try curl first, then wget, then R
        if command -v curl &> /dev/null; then
            curl -L https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R -o {output.visfuns} || exit 1
        elif command -v wget &> /dev/null; then
            wget -O {output.visfuns} https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R || exit 1
        else
            # Use R to download if curl/wget not available
            Rscript -e "download.file('https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R', '{output.visfuns}', method='libcurl')" || exit 1
        fi
        """

# Rule to install evaladmix from source (if not already installed)
rule install_evaladmix:
    output:
        touch(".snakemake/evaladmix_installed")
    log:
        "logs/evaladmix_install.log"
    conda:
        "../envs/evaladmix.yaml"
    shell:
        """
        set -e  # Exit on error
        set -x  # Print commands for debugging
        
        # Check if evalAdmix is already installed
        if ! command -v evalAdmix &> /dev/null; then
            echo "evalAdmix not found. Installing from source..."
            echo "Repository: https://github.com/GenisGE/evalAdmix"
            
            # Clone and build evaladmix from GitHub
            TMPDIR=$(mktemp -d)
            echo "Using temporary directory: $TMPDIR"
            cd $TMPDIR || {{ echo "Failed to cd to $TMPDIR"; exit 1; }}
            
            echo "Cloning evalAdmix repository from https://github.com/GenisGE/evalAdmix..."
            git clone https://github.com/GenisGE/evalAdmix.git || {{ echo "git clone failed"; exit 1; }}
            
            cd evalAdmix || {{ echo "Failed to cd to evalAdmix directory"; exit 1; }}
            
            echo "Building evalAdmix..."
            make || {{ echo "make failed. Checking Makefile..."; cat Makefile; exit 1; }}
            
            # Check if binary was created
            if [ ! -f evalAdmix ]; then
                echo "ERROR: evalAdmix binary not found after make"
                ls -la
                exit 1
            fi
            
            echo "Installing evalAdmix to $CONDA_PREFIX/bin"
            mkdir -p $CONDA_PREFIX/bin || {{ echo "Failed to create bin directory"; exit 1; }}
            cp evalAdmix $CONDA_PREFIX/bin/ || {{ echo "Failed to copy evalAdmix"; exit 1; }}
            
            # Verify installation
            if [ ! -f $CONDA_PREFIX/bin/evalAdmix ]; then
                echo "ERROR: evalAdmix not found in $CONDA_PREFIX/bin after installation"
                exit 1
            fi
            
            chmod +x $CONDA_PREFIX/bin/evalAdmix
            cd - || true
            rm -rf $TMPDIR
            echo "evalAdmix installed successfully at $CONDA_PREFIX/bin/evalAdmix"
        else
            echo "evalAdmix already installed at $(command -v evalAdmix)"
        fi
        """

rule evaladmix_admixture:
    """
    Run evalAdmix for ADMIXTURE results.
    """
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam,
        qfile = rules.admixture.output.Q,
        pfile = rules.admixture.output.P,
        install = rules.install_evaladmix.output
    output:
        corres = "results/{project}/admixture/evaladmix/{project}.admixture.K{k}.corres"
    log:
        "logs/{project}/evaladmix_admixture.K{k}.log"
    benchmark:
        "benchmarks/{project}/evaladmix_admixture.K{k}.txt"
    params:
        plink_prefix = "results/{project}/filtered_data/{project}.biallelic_snps",
        output_prefix = lambda wildcards, output: str(output.corres)
    conda:
        "../envs/evaladmix.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        evalAdmix \
            -plink {params.plink_prefix} \
            -fname {input.pfile} \
            -qname {input.qfile} \
            -o {params.output_prefix} \
            -P {threads} \
            2>&1 | tee {log}
        """

rule evaladmix_faststructure:
    """
    Run evalAdmix for fastStructure results.
    """
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam,
        qfile = rules.faststructure.output.meanQ,
        pfile = rules.faststructure.output.meanP,
        install = rules.install_evaladmix.output
    output:
        corres = "results/{project}/faststructure/evaladmix/{project}.faststructure.K{k}.corres"
    log:
        "logs/{project}/evaladmix_faststructure.K{k}.log"
    benchmark:
        "benchmarks/{project}/evaladmix_faststructure.K{k}.txt"
    params:
        plink_prefix = "results/{project}/filtered_data/{project}.biallelic_snps",
        output_prefix = lambda wildcards, output: str(output.corres)
    conda:
        "../envs/evaladmix.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        evalAdmix \
            -plink {params.plink_prefix} \
            -fname {input.pfile} \
            -qname {input.qfile} \
            -o {params.output_prefix} \
            -P {threads} \
            2>&1 | tee {log}
        """

rule structure_extract_pfile:
    """
    Extract P file (ancestral frequencies) from STRUCTURE output.
    This is a simplified extraction - for accurate results, a more sophisticated parser
    would extract actual allele frequencies from STRUCTURE output files.
    """
    input:
        structure_runs = lambda wildcards: expand(
            "results/{project}/structure/{project}.structure.K{k}.R{r}_f",
            project=wildcards.project,
            k=wildcards.k,
            r=range(1, config["projects"][wildcards.project]["parameters"]["structure"].get("replicates", 1) + 1)
        )
    output:
        pfile = "results/{project}/structure/evaladmix/{project}.structure.K{k}.P"
    log:
        "logs/{project}/structure_extract_pfile.K{k}.log"
    benchmark:
        "benchmarks/{project}/structure_extract_pfile.K{k}.txt"
    conda:
        "../envs/r-pophelper.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/extract_structure_pfile.R"

rule evaladmix_structure:
    """
    Run evalAdmix for STRUCTURE results.
    Note: This requires extracting P file from STRUCTURE output, which is complex.
    The current implementation uses a simplified approach.
    """
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam,
        qfile = rules.structure_summarize_pophelper.output.qmatrix,
        pfile = rules.structure_extract_pfile.output.pfile,
        install = rules.install_evaladmix.output
    output:
        corres = "results/{project}/structure/evaladmix/{project}.structure.K{k}.corres"
    log:
        "logs/{project}/evaladmix_structure.K{k}.log"
    benchmark:
        "benchmarks/{project}/evaladmix_structure.K{k}.txt"
    params:
        plink_prefix = "results/{project}/filtered_data/{project}.biallelic_snps",
        output_prefix = lambda wildcards, output: str(output.corres)
    conda:
        "../envs/evaladmix.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        set -e  # Exit on error
        mkdir -p $(dirname {output.corres})
        # Convert Q matrix from tab-delimited to space-delimited for evalAdmix
        tr '\t' ' ' < {input.qfile} > results/{wildcards.project}/structure/evaladmix/{wildcards.project}.structure.K{wildcards.k}.Q.tmp
        
        evalAdmix \
            -plink {params.plink_prefix} \
            -fname {input.pfile} \
            -qname results/{wildcards.project}/structure/evaladmix/{wildcards.project}.structure.K{wildcards.k}.Q.tmp \
            -o {params.output_prefix} \
            -P {threads} \
            2>&1 | tee {log}
        rm -f results/{wildcards.project}/structure/evaladmix/{wildcards.project}.structure.K{wildcards.k}.Q.tmp
        
        # Check if output file was created (evalAdmix appends .txt to the -o prefix)
        if [ ! -f {params.output_prefix}.txt ] && [ ! -f {output.corres} ]; then
            echo "ERROR: evalAdmix did not create output file" >&2
            echo "Expected: {output.corres} or {params.output_prefix}.txt" >&2
            ls -la $(dirname {output.corres})/ || true
            exit 1
        fi
        # If evalAdmix created .txt file, rename it to .corres
        if [ -f {params.output_prefix}.txt ] && [ ! -f {output.corres} ]; then
            mv {params.output_prefix}.txt {output.corres}
        fi
        """

rule plot_evaladmix_admixture:
    """
    Plot evalAdmix correlation matrix for ADMIXTURE.
    """
    input:
        corres = rules.evaladmix_admixture.output.corres,
        popmap = rules.generate_popmap.output,
        qfile = rules.admixture.output.Q,
        visfuns = rules.download_visfuns.output.visfuns
    output:
        pdf = "results/{project}/admixture/evaladmix/plots/{project}.admixture.K{k}.evaladmix.pdf",
        rds = "results/{project}/admixture/evaladmix/plots/{project}.admixture.K{k}.evaladmix.rds"
    log:
        "logs/{project}/plot_evaladmix_admixture.K{k}.log"
    benchmark:
        "benchmarks/{project}/plot_evaladmix_admixture.K{k}.txt"
    params:
        method = "ADMIXTURE",
        k = lambda wildcards: wildcards.k
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_evaladmix.R"

rule plot_evaladmix_faststructure:
    """
    Plot evalAdmix correlation matrix for fastStructure.
    """
    input:
        corres = rules.evaladmix_faststructure.output.corres,
        popmap = rules.generate_popmap.output,
        qfile = rules.faststructure.output.meanQ,
        visfuns = rules.download_visfuns.output.visfuns
    output:
        pdf = "results/{project}/faststructure/evaladmix/plots/{project}.faststructure.K{k}.evaladmix.pdf",
        rds = "results/{project}/faststructure/evaladmix/plots/{project}.faststructure.K{k}.evaladmix.rds"
    log:
        "logs/{project}/plot_evaladmix_faststructure.K{k}.log"
    benchmark:
        "benchmarks/{project}/plot_evaladmix_faststructure.K{k}.txt"
    params:
        method = "fastStructure",
        k = lambda wildcards: wildcards.k
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_evaladmix.R"

rule plot_evaladmix_structure:
    """
    Plot evalAdmix correlation matrix for STRUCTURE.
    """
    input:
        corres = rules.evaladmix_structure.output.corres,
        popmap = rules.generate_popmap.output,
        qfile = rules.structure_summarize_pophelper.output.qmatrix,
        visfuns = rules.download_visfuns.output.visfuns
    output:
        pdf = "results/{project}/structure/evaladmix/plots/{project}.structure.K{k}.evaladmix.pdf",
        rds = "results/{project}/structure/evaladmix/plots/{project}.structure.K{k}.evaladmix.rds"
    log:
        "logs/{project}/plot_evaladmix_structure.K{k}.log"
    benchmark:
        "benchmarks/{project}/plot_evaladmix_structure.K{k}.txt"
    params:
        method = "STRUCTURE",
        k = lambda wildcards: wildcards.k
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_evaladmix.R"


