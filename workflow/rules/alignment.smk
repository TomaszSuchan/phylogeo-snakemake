# Concatenated sequence alignment from ipyrad .loci, restricted to retained samples.
# Shared by IQ-TREE, trimAl+IQ-TREE, and rapidNJ.


rule loci_to_alignment:
    """
    Stream the ipyrad .loci file, keep only samples in samples_to_keep.txt, and
    write a concatenated relaxed-PHYLIP alignment plus FASTA.
    """
    input:
        loci=lambda wildcards: config["projects"][wildcards.project]["ipyrad_prefix"] + ".loci",
        samples=rules.filter_related_individuals.output.samples_to_keep,
    output:
        phy="results/{project}/filtered_data/{project}.loci.phy",
        fasta="results/{project}/filtered_data/{project}.loci.fasta",
        summary="results/{project}/filtered_data/{project}.loci_alignment_summary.tsv",
    log:
        "logs/{project}/loci_to_alignment.log",
    benchmark:
        "benchmarks/{project}/loci_to_alignment.txt",
    conda:
        "../envs/python.yaml",
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["loci_to_alignment"]["threads"],
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["loci_to_alignment"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["loci_to_alignment"]["runtime"],
    script:
        "../scripts/loci_to_alignment.py"
