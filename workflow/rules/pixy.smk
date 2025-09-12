rule prepare_invariant_vcf:
    input:
        loci = config["ipyrad_prefix"] + ".loci",
        vcf = rules.sort_vcf.output.vcf,
        index = rules.index_vcf.output.index
    output:
        invariant_vcf = config["analysis_name"] + "/filtered_data/invariant_sites.vcf.gz"
    log:
        config["analysis_name"] + "/logs/prepare_invariant_vcf.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/prepare_invariant_vcf.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        time = config["resources"]["default"]["runtime"]
    shell:
        """
        python workflow/scripts/extract_invariant_vcf.py {input.loci} -o {config['analysis_name']}/filtered_data/invariant_sites_only.vcf &> {log}
        bgzip {config['analysis_name']}/filtered_data/invariant_sites_only.vcf  >> {log} 2>&1
        bcftools index {config['analysis_name']}/filtered_data/invariant_sites_only.vcf.gz  >> {log} 2>&1
        bcftools concat {config['analysis_name']}/filtered_data/invariant_sites_only.vcf.gz {input.vcf} -Oz -a -o {output.invariant_vcf}  >> {log} 2>&1
        bcftools index {output.invariant_vcf}
        """

rule pixy:
    input:
        vcf = rules.prepare_invariant_vcf.output.invariant_vcf,
        popmap = "data/popmap.txt"
    output:
        pi = config["analysis_name"] + "/pixy/pi.txt",
        fst = config["analysis_name"] + "/pixy/fst.txt",
        dxy = config["analysis_name"] + "/pixy/dxy.txt"
    log:
        config["analysis_name"] + "/logs/pixy.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/pixy.txt"
    params:
        window_size = config["pixy"].get("window_size", 10000),
        output_folder = config["analysis_name"] + "/pixy/",
        output_prefix = ""
    conda:
        "../envs/pixy.yaml"
    threads: config["resources"]["pixy"]["threads"]
    resources:
        mem_mb = config["resources"]["pixy"]["mem_mb"],
        time = config["resources"]["pixy"]["runtime"]
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats pi fst dxy --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder} &> {log}
        """