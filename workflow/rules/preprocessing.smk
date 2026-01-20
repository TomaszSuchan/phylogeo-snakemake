import os

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
        samples_file="results/{project}/filtered_data/{project}.samples_subset.txt"
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
        vcf=temporary("results/{project}/filtered_data/{project}.vcf.gz")
    log:
        "logs/{project}/subset_vcf.log"
    benchmark:
        "benchmarks/{project}/subset_vcf.txt"
    params:
        samples=lambda wildcards: get_samples(wildcards),
        has_samples=lambda wildcards: "1" if get_samples(wildcards) else "0",
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
        if [ "{params.has_samples}" = "1" ]; then
            if [ "{params.skip_missing_filter}" = "1" ]; then
                # Skip F_MISSING filtering when f_missing >= 1
                bcftools view -S {input.samples_file} -Ou {input.vcf} | \
                bcftools +fill-tags - -Oz -o {output.vcf} -- -t NS,AC,AN,AF >> {log} 2>&1 || exit 1
            else
                # Apply F_MISSING filtering when f_missing < 1
                bcftools view -S {input.samples_file} -Ou {input.vcf} | \
                bcftools filter -e 'F_MISSING > {params.missing_threshold}' -Ou - | \
                bcftools +fill-tags - -Oz -o {output.vcf} -- -t NS,AC,AN,AF >> {log} 2>&1 || exit 1
            fi
        else
            if [ "{params.skip_missing_filter}" = "1" ]; then
                # Skip F_MISSING filtering when f_missing >= 1
                bcftools +fill-tags {input.vcf} -Oz -o {output.vcf} -- -t NS,AC,AN,AF >> {log} 2>&1 || exit 1
            else
                # Apply F_MISSING filtering when f_missing < 1
                bcftools filter -e 'F_MISSING > {params.missing_threshold}' -Ou {input.vcf} | \
                bcftools +fill-tags - -Oz -o {output.vcf} -- -t NS,AC,AN,AF >> {log} 2>&1 || exit 1
            fi
        fi
        """

# Rule to index subset VCF
rule index_subset_vcf:
    input:
        vcf=rules.subset_vcf.output.vcf
    output:
        index=temporary("results/{project}/filtered_data/{project}.vcf.gz.csi")
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

# Rule to filter related individuals using plink2 KING method
rule filter_related_individuals:
    input:
        vcf=rules.subset_vcf.output.vcf,
        index=rules.index_subset_vcf.output.index
    output:
        samples_to_keep="results/{project}/filtered_data/{project}.samples_to_keep.txt",
        king_table="results/{project}/filtered_data/{project}.king_table.tsv"
    log:
        "logs/{project}/filter_related_individuals.log"
    benchmark:
        "benchmarks/{project}/filter_related_individuals.txt"
    params:
        king_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("king_threshold", 0.0884),
        mac_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("mac_threshold", 1),
        plink_prefix="results/{project}/filtered_data/{project}.king_temp",
        enabled=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("enabled", False)
    conda:
        "../envs/plink.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("threads", 4)
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("mem_mb", 16000),
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("relatedness_filtering", {}).get("runtime", 120)
    shell:
        """
        if [ "{params.enabled}" = "True" ] || [ "{params.enabled}" = "true" ]; then
            # Convert VCF to PLINK format (pgen format supports multiallelic variants)
            plink2 --vcf {input.vcf} \
                   --make-pgen \
                   --out {params.plink_prefix} \
                   --allow-extra-chr 0 \
                   --double-id \
                   --mac {params.mac_threshold} >> {log} 2>&1 || exit 1
            
            # Generate full KING table (outputs .kin0 file)
            plink2 --pfile {params.plink_prefix} \
                   --make-king-table \
                   --out {params.plink_prefix} >> {log} 2>&1 || exit 1
            
            # Copy KING table to output location
            # plink2 outputs .kin0 file (or .kin if using rel-check mode)
            if [ -f {params.plink_prefix}.kin0 ]; then
                cp {params.plink_prefix}.kin0 {output.king_table} || exit 1
            elif [ -f {params.plink_prefix}.kin ]; then
                cp {params.plink_prefix}.kin {output.king_table} || exit 1
            else
                echo "Error: KING table file not found" >> {log}
                exit 1
            fi
            
            # Filter related individuals using --king-cutoff (greedy algorithm)
            plink2 --pfile {params.plink_prefix} \
                   --king-cutoff {params.king_threshold} \
                   --make-pgen \
                   --out {params.plink_prefix}.filtered >> {log} 2>&1 || exit 1
            
            # Extract samples to keep from filtered .psam file
            awk '{{print $2}}' {params.plink_prefix}.filtered.psam > {output.samples_to_keep} || exit 1
            
            # Clean up temporary files
            rm -f {params.plink_prefix}.pgen {params.plink_prefix}.pvar {params.plink_prefix}.psam
            rm -f {params.plink_prefix}.filtered.pgen {params.plink_prefix}.filtered.pvar {params.plink_prefix}.filtered.psam
            rm -f {params.plink_prefix}.log {params.plink_prefix}.filtered.log
        else
            # If filtering is disabled, create samples_to_keep with all samples from VCF
            case "{input.vcf}" in
                *.gz) bcftools query -l {input.vcf} > {output.samples_to_keep} 2>> {log} || exit 1 ;;
                *) bcftools query -l {input.vcf} > {output.samples_to_keep} 2>> {log} || exit 1 ;;
            esac
            # Create empty KING table
            touch {output.king_table}
        fi
        """

# Rule to update samples file based on relatedness filtering
rule update_samples_file:
    input:
        original_samples=rules.create_samples_file.output.samples_file,
        samples_to_keep=rules.filter_related_individuals.output.samples_to_keep
    output:
        filtered_samples="results/{project}/filtered_data/{project}.samples_subset_relatedness_filtered.txt"
    log:
        "logs/{project}/update_samples_file.log"
    params:
        enabled=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("enabled", False)
    shell:
        """
        if [ "{params.enabled}" = "True" ] || [ "{params.enabled}" = "true" ]; then
            # Filter original samples file to only include samples in samples_to_keep
            grep -F -f {input.samples_to_keep} {input.original_samples} > {output.filtered_samples} 2> {log} || exit 1
        else
            # If filtering is disabled, copy original samples file
            cp {input.original_samples} {output.filtered_samples} 2> {log} || exit 1
        fi
        """

# Rule to categorize removed individuals by relationship type
rule categorize_removed_individuals:
    input:
        king_table=rules.filter_related_individuals.output.king_table,
        samples_to_keep=rules.filter_related_individuals.output.samples_to_keep
    output:
        categorized="results/{project}/filtered_data/{project}.removed_individuals.tsv"
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
        # Don't reference rules.generate_popdata here because popdata.smk is included
        # after preprocessing.smk in workflow/Snakefile; use the file path so Snakemake
        # can resolve the producing rule later.
        popdata="results/{project}/indpopdata.txt"
    output:
        pdf="results/{project}/filtered_data/plots/{project}.removed_individuals_by_{group_by}.pdf",
        rds="results/{project}/filtered_data/plots/{project}.removed_individuals_by_{group_by}.rds"
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
        vcf=rules.subset_vcf.output.vcf,
        index=rules.index_subset_vcf.output.index,
        samples_file=rules.update_samples_file.output.filtered_samples
    output:
        vcf="results/{project}/filtered_data/{project}.vcf_filtered_relatedness.gz"
    log:
        "logs/{project}/subset_vcf_after_relatedness.log"
    benchmark:
        "benchmarks/{project}/subset_vcf_after_relatedness.txt"
    params:
        enabled=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("enabled", False)
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        if [ "{params.enabled}" = "True" ] || [ "{params.enabled}" = "true" ]; then
            # Subset VCF to only include filtered samples
            bcftools view -S {input.samples_file} -Ou {input.vcf} | \
            bcftools +fill-tags - -Oz -o {output.vcf} -- -t NS,AC,AN,AF >> {log} 2>&1 || exit 1
        else
            # If filtering is disabled, just copy the VCF
            cp {input.vcf} {output.vcf} 2>> {log} || exit 1
        fi
        """

# Rule to index VCF after relatedness filtering
rule index_vcf_after_relatedness:
    input:
        vcf=rules.subset_vcf_after_relatedness.output.vcf
    output:
        index="results/{project}/filtered_data/{project}.vcf_filtered_relatedness.gz.csi"
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

# Rule to select only biallelic SNPs with MAC>1 from VCF
rule select_biallelic_snps:
    input:
        vcf = lambda wildcards: (
            rules.subset_vcf_after_relatedness.output.vcf
            if config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("enabled", False)
            else rules.subset_vcf.output.vcf
        ),
        index = lambda wildcards: (
            rules.index_vcf_after_relatedness.output.index
            if config["projects"][wildcards.project]["parameters"].get("relatedness_filtering", {}).get("enabled", False)
            else rules.index_subset_vcf.output.index
        )
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

# Rule to thin VCF using the custom script
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

# Rule to convert thinned VCF to PLINK format
# renames chromosomes to "0" for downstream compatibility with "--allow-extra-chr 0"
rule vcf_to_plink:
    input:
        vcf = rules.thin_vcf.output.vcf
    output:
        bed = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.bed",
        bim = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.bim",
        fam = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.fam"
    log:
        "logs/{project}/vcf_to_plink.log"
    benchmark:
        "benchmarks/{project}/vcf_to_plink.txt"
    params:
        output_prefix = "results/{project}/filtered_data/{project}.biallelic_snps_thinned"
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

# Rule to convert thinned VCF to STRUCTURE format
rule vcf_to_structure:
    input:
        vcf = rules.thin_vcf.output.vcf
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

# Rule to filter VCF for missing data threshold
rule filter_missing_vcf:
    input:
        vcf = rules.thin_vcf.output.vcf
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
