import os

# Helper function to get thinning strategy for a project
def get_thinning_strategy(wildcards):
    """Get thinning strategy for a project, default to 'thinning'."""
    return config["projects"][wildcards.project]["parameters"].get("thinning_strategy", "thinning")

# Helper function to get output filename suffix based on thinning strategy
def get_thinning_output_suffix(wildcards):
    """Get output filename suffix based on thinning strategy."""
    strategy = get_thinning_strategy(wildcards)
    if strategy == "thinning":
        return "thinned"
    elif strategy == "ld_pruning":
        return "ld_pruned"
    elif strategy == "none":
        return "all_snps"
    else:
        raise ValueError(f"Unknown thinning strategy: {strategy}")

# Helper function to get the filtered VCF output file based on thinning strategy
def get_filtered_vcf_output(wildcards):
    """Get the filtered VCF output file based on thinning strategy."""
    suffix = get_thinning_output_suffix(wildcards)
    return f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_{suffix}.vcf.gz"

# Helper function to get PLINK bfile prefix based on thinning strategy
def get_plink_bfile_prefix(wildcards):
    """Get the PLINK bfile prefix based on thinning strategy."""
    suffix = get_thinning_output_suffix(wildcards)
    return f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_{suffix}"

# Helper function to get conda environment based on thinning strategy
def get_thinning_conda_env(wildcards):
    """Get conda environment based on thinning strategy."""
    strategy = get_thinning_strategy(wildcards)
    if strategy == "thinning":
        return "../envs/python.yaml"
    elif strategy == "ld_pruning":
        return "../envs/plink.yaml"
    else:  # none
        return "../envs/bcftools.yaml"  # for copying/compressing VCF

# Helper function to get samples list for a project
def get_samples(wildcards):
    """Get samples list for a project from comma-separated string, return empty list if not specified."""
    samples = config["projects"][wildcards.project].get("samples", "")
    
    # If samples is a string (comma-separated), split and convert to list
    if isinstance(samples, str):
        # Split by comma and strip whitespace, filter out empty strings
        return [s.strip() for s in samples.split(',') if s.strip()]
    
    # Return empty list if not a string
    return []

# Rule to sort input vcf
rule sort_vcf:
    input:
      vcf=lambda wildcards: (
          f"{config['projects'][wildcards.project]['ipyrad_prefix']}.vcf.gz"
          if os.path.exists(f"{config['projects'][wildcards.project]['ipyrad_prefix']}.vcf.gz")
          else f"{config['projects'][wildcards.project]['ipyrad_prefix']}.vcf"
      )
    output:
        vcf=temporary("results/{project}/filtered_data/{project}.raw_sorted.vcf.gz")
    log:
        "logs/{project}/sort_vcf.log"
    benchmark:
        "benchmarks/{project}/sort_vcf.txt"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcf-sort {input.vcf} 2> {log} | bgzip -c > {output.vcf}
        """

rule index_vcf:
    input:
        # Snakemake selects the first existing file from this list
        vcf=rules.sort_vcf.output.vcf
    output:
        index=temporary("results/{project}/filtered_data/{project}.raw_sorted.vcf.gz.csi")
    log:
        "logs/{project}/index_vcf.log"
    benchmark:
        "benchmarks/{project}/index_vcf.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools index -f {input.vcf} &> {log}
        """

# Rule to create samples text file
rule create_samples_file:
    input:
        vcf=rules.sort_vcf.output.vcf
    output:
        samples_file=temporary("results/{project}/filtered_data/{project}.samples_subset.txt")
    log:
        "logs/{project}/create_samples_file.log"
    benchmark:
        "benchmarks/{project}/create_samples_file.txt"
    params:
        samples=lambda wildcards: get_samples(wildcards),
        has_samples=lambda wildcards: "1" if get_samples(wildcards) else "0",
        samples_repr=lambda wildcards: repr(get_samples(wildcards))
    shell:
        """
        if [ "{params.has_samples}" = "1" ]; then
            # First, write the provided samples to a temp file
            python -c "samples = {params.samples_repr}; [print(s) for s in samples]" > {output.samples_file}.tmp
            # Extract samples from VCF in VCF order, then filter to match provided samples
            case "{input.vcf}" in
                *.gz) zgrep "^#CHROM" {input.vcf} | awk '{{for(i=10;i<=NF;i++) print $i}}' | grep -F -f {output.samples_file}.tmp > {output.samples_file} 2> {log} || (rm -f {output.samples_file}.tmp && exit 1) ;;
                *) grep "^#CHROM" {input.vcf} | awk '{{for(i=10;i<=NF;i++) print $i}}' | grep -F -f {output.samples_file}.tmp > {output.samples_file} 2> {log} || (rm -f {output.samples_file}.tmp && exit 1) ;;
            esac
            rm -f {output.samples_file}.tmp
        else
            case "{input.vcf}" in
                *.gz) zgrep "^#CHROM" {input.vcf} | awk '{{for(i=10;i<=NF;i++) print $i}}' > {output.samples_file} 2> {log} || exit 1 ;;
                *) grep "^#CHROM" {input.vcf} | awk '{{for(i=10;i<=NF;i++) print $i}}' > {output.samples_file} 2> {log} || exit 1 ;;
            esac
        fi
        """

# Rule to subset samples from VCF using bcftools
rule subset_vcf:
    input:
        vcf=rules.sort_vcf.output.vcf,
        index=rules.index_vcf.output.index,
        samples_file=rules.create_samples_file.output.samples_file
    output:
        # Temporary: only user-provided sample subset (no relatedness, no missingness filter)
        vcf=temporary("results/{project}/filtered_data/{project}.subset.vcf.gz")
    log:
        "logs/{project}/subset_vcf.log"
    benchmark:
        "benchmarks/{project}/subset_vcf.txt"
    params:
        samples=lambda wildcards: get_samples(wildcards),
        has_samples=lambda wildcards: "1" if get_samples(wildcards) else "0"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        if [ "{params.has_samples}" = "1" ]; then
            bcftools view -S {input.samples_file} -Ou {input.vcf} | \
            bcftools +fill-tags - -Oz -o {output.vcf} -- -t NS,AC,AN,AF >> {log} 2>&1 || exit 1
        else
            bcftools +fill-tags {input.vcf} -Oz -o {output.vcf} -- -t NS,AC,AN,AF >> {log} 2>&1 || exit 1
        fi
        """

# Rule to index subset VCF
rule index_subset_vcf:
    input:
        vcf=rules.subset_vcf.output.vcf
    output:
        index=temporary("results/{project}/filtered_data/{project}.subset.vcf.gz.csi")
    log:
        "logs/{project}/index_subset_vcf.log"
    benchmark:
        "benchmarks/{project}/index_subset_vcf.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools index -f {input.vcf} &> {log}
        """

# Rule to create samples list when relatedness filtering is disabled
rule create_samples_list_when_disabled:
    input:
        vcf=rules.subset_vcf.output.vcf,
        index=rules.index_subset_vcf.output.index
    output:
        samples_to_keep=temporary("results/{project}/filtered_data/{project}.samples_to_keep_disabled.txt")
    log:
        "logs/{project}/create_samples_list_when_disabled.log"
    benchmark:
        "benchmarks/{project}/create_samples_list_when_disabled.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        case "{input.vcf}" in
            *.gz) bcftools query -l {input.vcf} > {output.samples_to_keep} 2>> {log} || exit 1 ;;
            *) bcftools query -l {input.vcf} > {output.samples_to_keep} 2>> {log} || exit 1 ;;
        esac
        """

# Rule to convert VCF to PLINK format for KING analysis
rule convert_vcf_to_plink_for_king:
    input:
        vcf=rules.subset_vcf.output.vcf,
        index=rules.index_subset_vcf.output.index
    output:
        pgen=temporary("results/{project}/filtered_data/{project}.king_temp.pgen"),
        pvar=temporary("results/{project}/filtered_data/{project}.king_temp.pvar"),
        psam=temporary("results/{project}/filtered_data/{project}.king_temp.psam")
    log:
        "logs/{project}/convert_vcf_to_plink_for_king.log"
    benchmark:
        "benchmarks/{project}/convert_vcf_to_plink_for_king.txt"
    params:
        mac_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("mac_threshold", 1),
        plink_prefix="results/{project}/filtered_data/{project}.king_temp"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("threads", 4)
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("mem_mb", 16000),
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("runtime", 120)
    shell:
        """
        plink2 --vcf {input.vcf} \
               --make-pgen \
               --out {params.plink_prefix} \
               --allow-extra-chr 0 \
               --double-id \
               --mac {params.mac_threshold} >> {log} 2>&1 || exit 1
        """

# Rule to generate KING table (only runs when filtering is enabled)
rule generate_king_table:
    input:
        pgen=rules.convert_vcf_to_plink_for_king.output.pgen,
        pvar=rules.convert_vcf_to_plink_for_king.output.pvar,
        psam=rules.convert_vcf_to_plink_for_king.output.psam
    output:
        king_table=temporary("results/{project}/filtered_data/{project}.king_table_raw.tsv")
    log:
        "logs/{project}/generate_king_table.log"
    benchmark:
        "benchmarks/{project}/generate_king_table.txt"
    params:
        plink_prefix="results/{project}/filtered_data/{project}.king_temp"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("threads", 4)
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("mem_mb", 16000),
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("runtime", 120)
    shell:
        """
        plink2 --pfile {params.plink_prefix} \
               --make-king-table \
               --out {params.plink_prefix} >> {log} 2>&1 || exit 1
        
        # Copy KING table to output location
        if [ -f {params.plink_prefix}.kin0 ]; then
            cp {params.plink_prefix}.kin0 {output.king_table} || exit 1
        elif [ -f {params.plink_prefix}.kin ]; then
            cp {params.plink_prefix}.kin {output.king_table} || exit 1
        else
            echo "Error: KING table file not found" >> {log}
            exit 1
        fi

        # Cleanup PLINK side-effect files (they are not declared outputs)
        rm -f {params.plink_prefix}.kin0 {params.plink_prefix}.kin {params.plink_prefix}.log || true
        """

# Rule to filter related individuals using KING cutoff
rule filter_related_with_king:
    input:
        pgen=rules.convert_vcf_to_plink_for_king.output.pgen,
        pvar=rules.convert_vcf_to_plink_for_king.output.pvar,
        psam=rules.convert_vcf_to_plink_for_king.output.psam
    output:
        filtered_psam=temporary("results/{project}/filtered_data/{project}.king_temp.filtered.psam")
    log:
        "logs/{project}/filter_related_with_king.log"
    benchmark:
        "benchmarks/{project}/filter_related_with_king.txt"
    params:
        king_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("king_threshold", 0.0884),
        plink_prefix="results/{project}/filtered_data/{project}.king_temp"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("threads", 4)
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("mem_mb", 16000),
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("runtime", 120)
    shell:
        """
        plink2 --pfile {params.plink_prefix} \
               --king-cutoff {params.king_threshold} \
               --make-pgen \
               --out {params.plink_prefix}.filtered >> {log} 2>&1 || exit 1
        """

# Rule to extract samples to keep from filtered PLINK file
rule extract_samples_to_keep:
    input:
        filtered_psam=rules.filter_related_with_king.output.filtered_psam
    output:
        samples_to_keep=temporary("results/{project}/filtered_data/{project}.samples_to_keep_enabled.txt")
    log:
        "logs/{project}/extract_samples_to_keep.log"
    benchmark:
        "benchmarks/{project}/extract_samples_to_keep.txt"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    params:
        plink_prefix="results/{project}/filtered_data/{project}.king_temp"
    shell:
        """
        # plink2 .psam includes a header line; drop it and keep only IID (2nd column)
        awk 'NR>1 {{print $2}}' {input.filtered_psam} > {output.samples_to_keep} || exit 1

        # Cleanup PLINK side-effect files from --king-cutoff run
        rm -f {params.plink_prefix}.filtered.pgen {params.plink_prefix}.filtered.pvar {params.plink_prefix}.filtered.psam || true
        rm -f {params.plink_prefix}.filtered.log {params.plink_prefix}.log || true
        """

# Rule to combine outputs (samples_to_keep and king_table) - handles conditional logic
rule filter_related_individuals:
    input:
        samples_to_keep = lambda wildcards: (
            rules.extract_samples_to_keep.output.samples_to_keep
            if config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("enabled", False)
            else rules.create_samples_list_when_disabled.output.samples_to_keep
        ),
        king_table_raw=rules.generate_king_table.output.king_table
    output:
        samples_to_keep="results/{project}/filtered_data/{project}.samples_to_keep.txt",
        king_table="results/{project}/filtered_data/{project}.king_table.tsv"
    params:
        enabled=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("enabled", False)
    shell:
        """
        cp {input.samples_to_keep} {output.samples_to_keep}
        if [ "{params.enabled}" = "True" ] || [ "{params.enabled}" = "true" ]; then
            cp {input.king_table_raw} {output.king_table}
        else
            # Create empty KING table when filtering is disabled
            touch {output.king_table}
        fi
        """

#
# NOTE: We no longer write an extra "{project}.samples_subset_relatedness_filtered.txt" file.
# All downstream subsetting uses the canonical "{project}.samples_to_keep.txt".

# Rule to categorize removed individuals by relationship type
rule categorize_removed_individuals:
    input:
        king_table=rules.filter_related_individuals.output.king_table,
        samples_to_keep=rules.filter_related_individuals.output.samples_to_keep
    output:
        categorized="results/{project}/stats_samples/{project}.relatedness_filtered_samples.txt"
    log:
        "logs/{project}/categorize_removed_individuals.log"
    params:
        enabled=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("enabled", False),
        clone_threshold=0.354,  # duplicates/monozygotic twins
        first_degree_threshold=0.177,  # parents/siblings
        second_degree_threshold=0.0884  # half-siblings/grandparents
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        if [ "{params.enabled}" = "True" ] || [ "{params.enabled}" = "true" ]; then
            python workflow/scripts/categorize_removed_individuals.py \
                --king-table {input.king_table} \
                --samples-to-keep {input.samples_to_keep} \
                --output {output.categorized} \
                --clone-threshold {params.clone_threshold} \
                --first-degree-threshold {params.first_degree_threshold} \
                --second-degree-threshold {params.second_degree_threshold} &> {log} || exit 1
        else
            # If filtering is disabled, create empty file
            echo -e "individual\tmax_KING\tcategory" > {output.categorized} 2> {log} || exit 1
        fi
        """

# Rule to plot removed individuals barplot grouped by configurable column
# Only runs when relatedness filtering is enabled
rule plot_removed_individuals_barplot:
    input:
        removed_individuals=rules.categorize_removed_individuals.output.categorized,
        # Use indpopdata_all.txt which includes all individuals before filtering
        # (so removed individuals can be matched to their Site/Region)
        popdata="results/{project}/indpopdata_all.txt"
    output:
        pdf="results/{project}/stats_samples/plots/{project}.removed_individuals_by_{group_by}.pdf",
        rds="results/{project}/stats_samples/plots/{project}.removed_individuals_by_{group_by}.rds"
    log:
        "logs/{project}/plot_removed_individuals_barplot_{group_by}.log"
    wildcard_constraints:
        group_by="(?!sorted-)(?!grouped-).*"  # Exclude values starting with "sorted-" or "grouped-"
    params:
        enabled=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("enabled", False),
        group_by=lambda wildcards: wildcards.group_by if wildcards.group_by != "none" else None
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_removed_individuals_barplot.R"

# Rule to subset VCF after relatedness filtering (if enabled)
rule subset_vcf_after_relatedness:
    input:
        # Take the user-subset VCF, then apply relatedness-filtered sample list,
        # then apply missingness filtering (after relatedness), and write final VCF.
        vcf=rules.subset_vcf.output.vcf,
        index=rules.index_subset_vcf.output.index,
        samples_file=rules.filter_related_individuals.output.samples_to_keep
    output:
        vcf="results/{project}/filtered_data/{project}.filtered.vcf.gz"
    log:
        "logs/{project}/subset_vcf_after_relatedness.log"
    benchmark:
        "benchmarks/{project}/subset_vcf_after_relatedness.txt"
    params:
        missing_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("vcf_filtering", {}).get("f_missing", 1),
        skip_missing_filter=lambda wildcards: "1" if config["projects"][wildcards.project]["parameters"].get("vcf_filtering", {}).get("f_missing", 1) >= 1 else "0"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        if [ "{params.skip_missing_filter}" = "1" ]; then
            bcftools view -S {input.samples_file} -Ou {input.vcf} | \
            bcftools +fill-tags - -Oz -o {output.vcf} -- -t NS,AC,AN,AF >> {log} 2>&1 || exit 1
        else
            bcftools view -S {input.samples_file} -Ou {input.vcf} | \
            bcftools filter -e 'F_MISSING > {params.missing_threshold}' -Ou - | \
            bcftools +fill-tags - -Oz -o {output.vcf} -- -t NS,AC,AN,AF >> {log} 2>&1 || exit 1
        fi
        """

# Rule to index VCF after relatedness filtering
rule index_vcf_after_relatedness:
    input:
        vcf=rules.subset_vcf_after_relatedness.output.vcf
    output:
        index="results/{project}/filtered_data/{project}.filtered.vcf.gz.csi"
    log:
        "logs/{project}/index_vcf_after_relatedness.log"
    benchmark:
        "benchmarks/{project}/index_vcf_after_relatedness.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools index -f {input.vcf} &> {log}
        """

## Rules to filter VCF for biallelic unlinked markers
# Rule to select only biallelic SNPs with MAC>1 from VCF
rule select_biallelic_snps:
    input:
        vcf=rules.subset_vcf_after_relatedness.output.vcf,
        index=rules.index_vcf_after_relatedness.output.index
    output:
        biallelic_vcf = "results/{project}/filtered_data/{project}.biallelic_snps.vcf.gz"
    log:
        "logs/{project}/select_biallelic_snps.log"
    benchmark:
        "benchmarks/{project}/select_biallelic_snps.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools view -i 'MAC > 1' -m2 -M2 \
        -v snps {input.vcf} \
        -Oz -o {output.biallelic_vcf} &> {log}
        """

# Rule to thin VCF using Python script (strategy: "thinning")
rule thin_vcf:
    input:
        vcf = rules.select_biallelic_snps.output.biallelic_vcf
    output:
        vcf = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.vcf.gz"
    log:
        "logs/{project}/thin_vcf.log"
    benchmark:
        "benchmarks/{project}/thin_vcf.txt"
    params:
        min_coverage = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("min_coverage", 0),
        method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("method", "max_coverage"),
        ties = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("ties", "random"),
        ns_tag = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("ns_tag", "NS"),
        id_pattern = lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf_thinning"].get("id_pattern", r"loc(\d+)_")
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        python workflow/scripts/thin_ipyrad_vcf.py \
            --vcf {input.vcf} \
            --out {output.vcf} \
            --min-coverage {params.min_coverage} \
            --method {params.method} \
            --ties {params.ties} \
            --ns-tag {params.ns_tag} \
            --id-pattern '{params.id_pattern}' &> {log}
        """

# LD pruning rules are now in ld_pruning.smk
# The final output rule is ld_prune_compress

# Rule to pass through all SNPs without filtering (strategy: "none")
rule pass_through_vcf:
    input:
        vcf = rules.select_biallelic_snps.output.biallelic_vcf
    output:
        vcf = "results/{project}/filtered_data/{project}.biallelic_snps_all_snps.vcf.gz"
    log:
        "logs/{project}/pass_through_vcf.log"
    benchmark:
        "benchmarks/{project}/pass_through_vcf.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default-long"]["runtime"]
    shell:
        """
        # Just copy the biallelic VCF (it's already compressed)
        cp {input.vcf} {output.vcf} || exit 1
        """

## Rules for exporting data from VCF to other formats
# Rule to convert filtered VCF to PLINK format
# renames chromosomes to "0" for downstream compatibility with "--allow-extra-chr 0"
rule vcf_to_plink:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        bed = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.bed",
        bim = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.bim",
        fam = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.fam"
    log:
        "logs/{project}/vcf_to_plink.log"
    benchmark:
        "benchmarks/{project}/vcf_to_plink.txt"
    params:
        output_prefix = lambda wildcards: f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_{get_thinning_output_suffix(wildcards)}"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --out {params.output_prefix} \
              --allow-extra-chr 0 \
              --double-id &> {log}
        # Move files to fixed output names (always thinned for compatibility)
        mv {params.output_prefix}.bed {output.bed}
        mv {params.output_prefix}.bim {output.bim}
        mv {params.output_prefix}.fam {output.fam}
        """

# Rule to convert filtered VCF to STRUCTURE format
rule vcf_to_structure:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        str = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.str"
    log:
        "logs/{project}/vcf_to_structure.log"
    benchmark:
        "benchmarks/{project}/vcf_to_structure.txt"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/vcf_to_structure.py \
            --vcf {input.vcf} \
            --out {output.str} &> {log}
        """

## Rules used for genertaing VCFs with different data thresholds - for checking PCA robustness to missing data
# Rule to filter VCF for missing data threshold
rule filter_missing_vcf:
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        vcf = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.vcf.gz"
    log:
        "logs/{project}/filter_missing_vcf_{miss}.log"
    benchmark:
        "benchmarks/{project}/filter_missing_vcf_{miss}.txt"
    params:
        miss = lambda wildcards: wildcards.miss
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools view -i 'F_MISSING<{params.miss}' {input.vcf} -Oz -o {output.vcf} &> {log}
        """

# Rule to convert to PLINK VCFs filtered for missing data threshold
# renames chromosomes to "0" for downstream compatibility with "--allow-extra-chr 0"
rule missing_vcf_to_plink:
    input:
        vcf = rules.filter_missing_vcf.output.vcf
    output:
        bed = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.bed",
        bim = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.bim",
        fam = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.fam"
    log:
        "logs/{project}/missing_vcf_to_plink_{miss}.log"
    benchmark:
        "benchmarks/{project}/missing_vcf_to_plink_{miss}.txt"
    params:
        miss = lambda wildcards: wildcards.miss,
        output_prefix = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}"
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --out {params.output_prefix} \
              --allow-extra-chr 0 \
              --double-id &> {log}
        """
