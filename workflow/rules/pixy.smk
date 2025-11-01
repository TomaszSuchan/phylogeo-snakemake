rule prepare_invariant_vcf:
    input:
        loci = lambda wildcards: config["projects"][wildcards.project]["ipyrad_prefix"] + ".loci",
        index = rules.index_vcf.output.index
    output:
        invariant_vcf = "{project}/filtered_data/invariant_sites.vcf"
    log:
        "{project}/logs/prepare_invariant_vcf.log"
    benchmark:
        "{project}/benchmarks/prepare_invariant_vcf.txt"
    conda:
        "../envs/vcfpy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/extract_invariant_vcf.py {input.loci} -o {output.invariant_vcf} &> {log}
        """

rule prepare_invariant_vcf_gz:
    input:
        vcf = rules.prepare_invariant_vcf.output.invariant_vcf
    output:
        vcf = "{project}/filtered_data/invariant_sites.vcf.gz"
    log:
        "{project}/logs/prepare_invariant_vcf_gz.log"
    benchmark:
        "{project}/benchmarks/prepare_invariant_vcf_gz.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bgzip {input.vcf} &> {log}
        """

rule prepare_invariant_vcf_gz_index:
    input:
        vcf = rules.prepare_invariant_vcf_gz.output.vcf
    output:
        index = "{project}/filtered_data/invariant_sites.vcf.gz.csi"
    log:
        "{project}/logs/prepare_invariant_vcf_gz_index.log"
    benchmark:
        "{project}/benchmarks/prepare_invariant_vcf_gz_index.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools index {input.vcf}  &> {log}
        """

rule merge_invariant_sites:
    input:
        vcf_var = rules.sort_vcf.output.vcf,
        vcf_var_index = rules.index_vcf.output.index,
        vcf_inv = rules.prepare_invariant_vcf_gz.output.vcf,
        vcf_inv_index = rules.prepare_invariant_vcf_gz_index.output.index
    output:
        merged_vcf = "{project}/filtered_data/merged_invariant_sites.vcf.gz",
        index = "{project}/filtered_data/merged_invariant_sites.vcf.gz.csi"
    log:
        "{project}/logs/merge_invariant_sites.log"
    benchmark:
        "{project}/benchmarks/merge_invariant_sites.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        bcftools concat {input.vcf_inv} {input.vcf_var} -Oz -a -o {output.merged_vcf}  >> {log} 2>&1
        bcftools index {output.merged_vcf}
        """

rule pixy:
    input:
        vcf = rules.merge_invariant_sites.output.merged_vcf,
        popmap = rules.generate_popmap.output
    output:
        pi = "{project}/pixy/pixy_pi.txt",
        fst = "{project}/pixy/pixy_fst.txt",
        dxy = "{project}/pixy/pixy_dxy.txt"
    log:
        "{project}/logs/pixy.log"
    benchmark:
        "{project}/benchmarks/pixy.txt"
    params:
        window_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["pixy"].get("window_size", 10000),
        output_folder = "{project}/pixy/",
        output_prefix = ""
    conda:
        "../envs/pixy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pixy"]["runtime"]
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats pi fst dxy --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder} &> {log}
        """

rule pixy_summary:
    input:
        pi = rules.pixy.output.pi,
        fst = rules.pixy.output.fst,
        dxy = rules.pixy.output.dxy
    output:
        pi = "{project}/pixy/pixy_pi-summary.txt",
        fst = "{project}/pixy/pixy_fst-summary.txt",
        dxy = "{project}/pixy/pixy_dxy-summary.txt"
    log:
        "{project}/logs/pixy_summary.log"
    benchmark:
        "{project}/benchmarks/pixy_summary.txt"
    conda:
        "../envs/pandas.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/pixy_summary.py \
          --pi {input.pi} \
          --output-pi {output.pi} \
          --fst {input.fst} \
          --output-fst {output.fst} \
          --dxy {input.dxy} \
          --output-dxy {output.dxy} &> {log}
        """