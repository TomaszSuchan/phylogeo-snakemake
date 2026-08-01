"""
Folded site frequency spectra with easySFS (https://github.com/isaacovercast/easySFS).

Prep for Stairway Plot 2. Outputs live under
results/{project}/stairwayplot2/easysfs/{grouping}/.
Config: parameters.stairwayplot2.easysfs.group_by (or legacy population_column).
SNPs are taken from the all-sites VCF (no MAC/MAF filters from select_biallelic_snps)
so rare variants are retained and L is counted over the same samples and loci as the SNPs.

Moments SFS builds its own joint SFS under results/{project}/moments/sfs/ (see moments.smk).
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
        samples_dir=directory("results/{project}/stairwayplot2/easysfs/{grouping}/samples"),
        populations="results/{project}/stairwayplot2/easysfs/{grouping}/{project}.easysfs_populations.tsv",
    params:
        mode="column",
        population_column=lambda wildcards: wildcards.grouping,
        min_individuals=lambda wildcards: _stairwayplot2_easysfs_cfg(wildcards.project).get("min_individuals", 10),
        project_to_n_diploids=lambda wildcards: _stairwayplot2_easysfs_cfg(wildcards.project).get("project_to_n_diploids", None),
        max_project_diploids=lambda wildcards: _stairwayplot2_easysfs_cfg(wildcards.project).get("max_project_diploids", 40),
    log:
        "logs/{project}/easysfs_prepare_samples.{grouping}.log"
    benchmark:
        "benchmarks/{project}/easysfs_prepare_samples_{grouping}.txt"
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
        vcf="results/{project}/stairwayplot2/easysfs/{grouping}/vcf/{project}.{stratum}.biallelic.vcf.gz",
        index="results/{project}/stairwayplot2/easysfs/{grouping}/vcf/{project}.{stratum}.biallelic.vcf.gz.csi",
    params:
        samples_file="results/{project}/stairwayplot2/easysfs/{grouping}/samples/{project}.{stratum}.samples.txt",
    log:
        "logs/{project}/easysfs_subset_vcf.{grouping}.{stratum}.log"
    benchmark:
        "benchmarks/{project}/easysfs_subset_vcf_{grouping}_{stratum}.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("stairwayplot2", {}).get("mem_mb", 8000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("stairwayplot2", {}).get("runtime", 120),
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
        L="results/{project}/stairwayplot2/easysfs/{grouping}/{project}.{stratum}.L.txt",
    params:
        samples_file="results/{project}/stairwayplot2/easysfs/{grouping}/samples/{project}.{stratum}.samples.txt",
    log:
        "logs/{project}/easysfs_count_L.{grouping}.{stratum}.log"
    benchmark:
        "benchmarks/{project}/easysfs_count_L_{grouping}_{stratum}.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("stairwayplot2", {}).get("mem_mb", 8000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("stairwayplot2", {}).get("runtime", 120),
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
        sfs_dir=directory("results/{project}/stairwayplot2/easysfs/{grouping}/sfs/{project}.{stratum}"),
    params:
        popmap="results/{project}/stairwayplot2/easysfs/{grouping}/samples/{project}.{stratum}.popmap.txt",
        proj_file="results/{project}/stairwayplot2/easysfs/{grouping}/samples/{project}.{stratum}.proj.txt",
    log:
        "logs/{project}/easysfs_run.{grouping}.{stratum}.log"
    benchmark:
        "benchmarks/{project}/easysfs_run_{grouping}_{stratum}.txt"
    conda:
        "../envs/easysfs.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("stairwayplot2", {}).get("mem_mb", 32000),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("stairwayplot2", {}).get("runtime", 120),
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
