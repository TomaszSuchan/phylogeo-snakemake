rule generate_popmap:
    input:
        vcf=rules.sort_vcf.output.vcf
    output:
        popmap="results/{project}/popmap.txt"
    params:
        popmap=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popmap", ""),
        separator=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popseparator", "-")
    log:
        "logs/{project}/generate_popmap.log"
    benchmark:
        "benchmarks/{project}/generate_popmap.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        r"""
        exec &> >(tee -a {log})

        if [ -n "{params.popmap}" ] && [ -f "{params.popmap}" ]; then
            echo "Checking popmap against VCF samples..."

            # Get list of samples from VCF
            bcftools query -l {input.vcf} > {output.popmap}.vcf_samples.tmp
            vcf_sample_count=$(wc -l < {output.popmap}.vcf_samples.tmp)
            echo "Found $vcf_sample_count samples in VCF"

            # Build lookup table from popmap (sample -> population)
            # Then iterate through VCF samples in order and output matching entries
            awk 'NR==FNR{{popmap[$1]=$2; next}} {{if ($1 in popmap) {{print $1 "\t" popmap[$1]}} else {{missing[$1]=1}}}} END{{for (s in missing) {{print s}}}}' \
                {params.popmap} {output.popmap}.vcf_samples.tmp > {output.popmap} 2> {output.popmap}.missing.tmp

            # Check if any VCF samples are missing from popmap
            if [ -s {output.popmap}.missing.tmp ]; then
                missing_count=$(wc -l < {output.popmap}.missing.tmp)
                echo ""
                echo "ERROR: $missing_count sample(s) found in VCF but not in the provided popmap:"
                echo "----------------------------------------"
                cat {output.popmap}.missing.tmp
                echo "----------------------------------------"
                echo ""
                echo "Please update the popmap file to include all VCF samples."
                rm {output.popmap}.vcf_samples.tmp {output.popmap}.missing.tmp
                exit 1
            fi

            matched_count=$(wc -l < {output.popmap})
            echo "Successfully matched $matched_count samples"
            echo "Popmap generated with samples in VCF order"

            rm {output.popmap}.vcf_samples.tmp {output.popmap}.missing.tmp
        else
            echo "No popmap provided, generating from VCF sample names..."
            bcftools query -l {input.vcf} | \
            awk -v sep="{params.separator}" '{{
                split($0,a,sep);
                pop = (length(a)>1) ? a[1] : "UNKNOWN";
                print $0 "\t" pop
            }}' > {output.popmap}
            sample_count=$(wc -l < {output.popmap})
            echo "Generated popmap for $sample_count samples"
        fi
        """

rule generate_popdata:
    input:
        popmap=rules.generate_popmap.output.popmap
    output:
        indpopdata="results/{project}/indpopdata.txt"
    params:
        popdata=lambda wildcards: config["projects"][wildcards.project]["parameters"]["popdata"],
    log:
        "logs/{project}/generate_popdata.log"
    benchmark:
        "benchmarks/{project}/generate_popdata.txt"
    conda:
        "../envs/pandas.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/generate_popdata.py"



rule calculate_missing_indv:
    input:
        vcf = rules.thin_vcf.output.vcf
    output:
        imiss = "results/{project}/filtered_data/{project}.biallelic_snps_thinned.imiss"
    log:
        "logs/{project}/calculate_missing_indv.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_thinned"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --missing-indv \
                 --out {params.out_prefix} &> {log}
        """

rule calculate_missing_indv_miss:
    input:
        vcf = rules.filter_missing_vcf.output.vcf
    output:
        imiss = "results/{project}/filtered_data/{project}.biallelic_snps_thinned_miss{miss}.imiss"
    log:
        "logs/{project}/calculate_missing_indv_miss_{miss}.log"
    params:
        out_prefix = lambda wildcards: f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_thinned_miss{wildcards.miss}"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --missing-indv \
                 --out {params.out_prefix} &> {log}
        """