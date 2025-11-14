import hashlib

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

def generate_unique_seed(wildcards):
    """Generate a deterministic but unique seed based on wildcards"""
    # Combine wildcards to create a unique identifier for this specific job
    # This is deterministic - same wildcards always produce same seed
    job_id = f"{wildcards.project}_{wildcards.k}_{wildcards.r}"
    # Create a hash and convert to a 10-digit integer
    hash_obj = hashlib.sha256(job_id.encode())
    seed = int(hash_obj.hexdigest()[:10], 16) % 10000000000
    return seed

# Rule for structure
rule structure:
    input:
        ustr = rules.vcf_to_structure.output.str
    output:
        stroutput = "results/{project}/structure/{project}.structure.K{k}.R{r}_f"
    benchmark:
        "benchmarks/{project}/structure.K{k}.R{r}.txt"
    params:
        mainparams = "data/mainparams",
        extraparams = "data/extraparams",
        nind = lambda wildcards, input: get_ustr_info(input.ustr)[0],
        nloci = lambda wildcards, input: get_ustr_info(input.ustr)[1],
        basename = lambda wildcards: f"results/{wildcards.project}/structure/{wildcards.project}.structure.K{wildcards.k}.R{wildcards.r}",
        seed = lambda wildcards: generate_unique_seed(wildcards)
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
            -o {params.basename} \
            -D {params.seed}
        """

# Summarize STRUCTURE runs with pophelper and output Q matrix for mapmixture
rule structure_summarize_pophelper:
    input:
        lambda wildcards: expand(
            "results/{project}/structure/{project}.structure.K{k}.R{r}_f",
            project=wildcards.project,
            k=wildcards.k,
            r=range(1, config["projects"][wildcards.project]["parameters"]["structure"].get("replicates", 1) + 1)
        )
    output:
        qmatrix = "results/{project}/structure/{project}.structure.K{k}.Qmatrix.txt",
        aligned_runs = "results/{project}/structure/{project}.structure.K{k}.aligned_runs.rds"
    log:
        "logs/{project}/structure_summarize_pophelper.K{k}.log"
    conda:
        "../envs/r-pophelper.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/structure_summarize_pophelper.R"

# Plot aligned K values
rule plot_aligned_k:
    input:
        lambda wildcards: expand(
            "results/{project}/structure/{project}.structure.K{k}.R{r}_f",
            project=wildcards.project,
            k=config["projects"][wildcards.project]["parameters"]["k_values"],
            r=range(1, config["projects"][wildcards.project]["parameters"]["structure"].get("replicates", 1) + 1)
        )
    output:
        "results/{project}/structure/{project}.K_aligned.pdf"
    log:
        "logs/{project}/plot_K_aligned.log"
    conda:
        "../envs/r-pophelper.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_k_aligned.R"