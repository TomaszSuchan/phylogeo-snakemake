# Function to parse .ustr file and get NIND and NLOCI
def get_ustr_info(ustr_file):
    """Parse .ustr file to get number of individuals and loci"""
    with open(ustr_file, 'r') as f:
        lines = f.readlines()
    
    nind = 0
    nloci = 0
    for i, line in enumerate(lines):
        if line.strip():  # skip empty lines
            if nloci == 0:  # first data line
                # Count loci from first line (excluding individual name and pop info)
                parts = line.strip().split('\t')
                nloci = len(parts) - 6  # subtract individual name and population columns
            nind += 1
    
    nind = nind // 2  # each individual has 2 lines in STRUCTURE format
    return nind, nloci


# Rule for structure
rule structure:
    input:
        ustr = rules.vcf_to_structure.output.str
    output:
        stroutput = "{project}/structure/structure.K{k}.R{r}_f"
    benchmark:
        "{project}/benchmarks/structure.K{k}.R{r}.txt"
    params:
        mainparams = "data/mainparams",
        extraparams = "data/extraparams",
        nind = lambda wildcards, input: get_ustr_info(input.ustr)[0],
        nloci = lambda wildcards, input: get_ustr_info(input.ustr)[1],
        basename = lambda wildcards: f"{config['projects'][wildcards.project]}/structure/structure.K{wildcards.k}.R{wildcards.r}"
    conda:
        "../envs/structure.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["structure"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["structure"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["structure"]["runtime"]
    shell:
        """
        structure \
            -K {wildcards.k} \
            -L {params.nloci} \
            -N {params.nind} \
            -m {params.mainparams} \
            -e {params.extraparams} \
            -i {input.ustr} \
            -o {params.basename}
        """

# Plot aligned K values
rule plot_aligned_k:
    input:
        lambda wildcards: expand(
            "{project}/structure/structure.K{k}.R{r}_f",
            project=wildcards.project,
            k=config["projects"][wildcards.project]["parameters"]["k_values"],
            r=range(1, config["projects"][wildcards.project]["parameters"]["structure"].get("replicates", 1) + 1)
        )
    output:
        "{project}/structure/plots/K{k}_aligned.pdf"
    log:
        "{project}/logs/plot_k_aligned_K{k}.log"
    conda:
        "../envs/r-pophelper.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_k_aligned.R"