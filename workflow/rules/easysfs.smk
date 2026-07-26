"""
Folded site frequency spectra with easySFS (https://github.com/isaacovercast/easySFS).

Shared prep for SFS-based methods (Stairway Plot 2, and dadi/moments/fastsimcoal2 if added
later). SNPs are taken from the all-sites VCF, i.e. without the MAC/MAF filters applied in
select_biallelic_snps, so that rare variants are retained and the number of callable sites (L)
is counted over exactly the same samples and loci as the SNPs.
"""


rule easysfs_install:
    output:
        easysfs=".snakemake/easySFS/easySFS.py",
    params:
        url="https://github.com/isaacovercast/easySFS.git",
        clone_dir=".snakemake/easySFS",
    log:
        "logs/easysfs_install.log"
    conda:
        "../envs/easysfs.yaml"
    localrule: True
    shell:
        r"""
        set -euo pipefail
        # git clone refuses a non-empty target, so drop any earlier clone first.
        rm -rf "{params.clone_dir}"
        git clone --depth 1 "{params.url}" "{params.clone_dir}" > {log} 2>&1
        """


rule easysfs_prepare_samples:
    """Per-population sample lists, easySFS popmaps and --proj values."""
    input:
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        samples_dir=directory("results/{project}/easysfs/samples"),
        populations="results/{project}/easysfs/{project}.easysfs_populations.tsv",
    params:
        population_column=lambda wildcards: config["projects"][wildcards.project]["parameters"]["easysfs"].get("population_column", "Site"),
        min_individuals=lambda wildcards: config["projects"][wildcards.project]["parameters"]["easysfs"].get("min_individuals", 10),
        project_to_n_diploids=lambda wildcards: config["projects"][wildcards.project]["parameters"]["easysfs"].get("project_to_n_diploids", None),
        max_project_diploids=lambda wildcards: config["projects"][wildcards.project]["parameters"]["easysfs"].get("max_project_diploids", 40),
    log:
        "logs/{project}/easysfs_prepare_samples.log"
    benchmark:
        "benchmarks/{project}/easysfs_prepare_samples.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/easysfs_prepare_samples.py"


rule easysfs_subset_vcf:
    """Per-population biallelic SNPs, unfiltered by allele frequency."""
    input:
        vcf=rules.prepare_invariant_vcf_gz.output.vcf,
        vcf_index=rules.prepare_invariant_vcf_gz_index.output.index,
        samples=rules.easysfs_prepare_samples.output.samples_dir,
    output:
        vcf="results/{project}/easysfs/vcf/{project}.{stratum}.biallelic.vcf.gz",
        index="results/{project}/easysfs/vcf/{project}.{stratum}.biallelic.vcf.gz.csi",
    params:
        samples_file="results/{project}/easysfs/samples/{project}.{stratum}.samples.txt",
    log:
        "logs/{project}/easysfs_subset_vcf.{stratum}.log"
    benchmark:
        "benchmarks/{project}/easysfs_subset_vcf_{stratum}.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("easysfs", {}).get("mem_mb", 8000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("easysfs", {}).get("runtime", 120),
    shell:
        r"""
        set -euo pipefail
        bcftools view \
            -S {params.samples_file} \
            -m2 -M2 -v snps \
            -Oz -o {output.vcf} \
            {input.vcf} > {log} 2>&1
        bcftools index -f {output.vcf} >> {log} 2>&1
        """


rule easysfs_count_L:
    """Total callable sites (variant + invariant) for easySFS --total-length."""
    input:
        vcf=rules.prepare_invariant_vcf_gz.output.vcf,
        vcf_index=rules.prepare_invariant_vcf_gz_index.output.index,
        samples=rules.easysfs_prepare_samples.output.samples_dir,
    output:
        L="results/{project}/easysfs/{project}.{stratum}.L.txt",
    params:
        samples_file="results/{project}/easysfs/samples/{project}.{stratum}.samples.txt",
    log:
        "logs/{project}/easysfs_count_L.{stratum}.log"
    benchmark:
        "benchmarks/{project}/easysfs_count_L_{stratum}.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("easysfs", {}).get("mem_mb", 8000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("easysfs", {}).get("runtime", 120),
    shell:
        r"""
        set -euo pipefail
        bcftools view -S {params.samples_file} -H {input.vcf} 2> {log} | wc -l > {output.L}
        """


rule easysfs_run:
    """Folded 1D SFS in dadi and fastsimcoal2 format, down-projected to --proj gene copies."""
    input:
        easysfs=rules.easysfs_install.output.easysfs,
        vcf=rules.easysfs_subset_vcf.output.vcf,
        L=rules.easysfs_count_L.output.L,
        samples=rules.easysfs_prepare_samples.output.samples_dir,
    output:
        sfs_dir=directory("results/{project}/easysfs/sfs/{project}.{stratum}"),
    params:
        popmap="results/{project}/easysfs/samples/{project}.{stratum}.popmap.txt",
        proj_file="results/{project}/easysfs/samples/{project}.{stratum}.proj.txt",
    log:
        "logs/{project}/easysfs_run.{stratum}.log"
    benchmark:
        "benchmarks/{project}/easysfs_run_{stratum}.txt"
    conda:
        "../envs/easysfs.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("easysfs", {}).get("mem_mb", 32000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("easysfs", {}).get("runtime", 120),
    shell:
        r"""
        set -euo pipefail
        python {input.easysfs} \
            -i {input.vcf} \
            -p {params.popmap} \
            -o {output.sfs_dir} \
            --proj "$(cat {params.proj_file})" \
            --total-length "$(cat {input.L})" \
            -a -f -y \
            > {log} 2>&1
        """
