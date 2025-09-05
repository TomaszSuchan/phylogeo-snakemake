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
        stroutput = "results/structure/structure.K{k}.R{r}_f"
    params:
        mainparams = "data/mainparams",
        extraparams = "data/extraparams",
        nind = lambda wildcards, input: get_ustr_info(input.ustr)[0],
        nloci = lambda wildcards, input: get_ustr_info(input.ustr)[1],
        basename = lambda wildcards: f"results/structure/structure.K{wildcards.k}.R{wildcards.r}"
    conda:
        "../envs/structure.yaml"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = "20:00:00"
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

# Rule to run structure for all K values  
rule run_structure:
    input:
        expand("results/structure/structure.K{k}.R{r}_f", k=config["k_values"], r=list(range(1, config["structure"].get("replicates", 1) + 1)))