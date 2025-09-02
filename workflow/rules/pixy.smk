rule prepare_invariant_vcf:
    input:
        loci = config["ipyrad_prefix"] + ".loci",
        vcf = rules.sort_vcf.output.vcf,
        index = rules.index_vcf.output.index
    output:
        invariant_vcf = "filtered_data/invariant_sites.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb = 4000,
        time = "10:00"
    shell:
        """
        python workflow/scripts/extract_invariant_vcf.py {input.loci} -o filtered_data/invariant_sites.vcf
        bgzip filtered_data/invariant_sites.vcf
        bcftools index invariant_sites.vcf.gz
        bcftools concat invariant_sites.vcf.gz {input.vcf} -Oz -a -o {output.invariant_vcf}
        bcftools index {output.invariant_vcf}
        """

rule pixy:
    input:
        vcf = rules.prepare_invariant_vcf.output.invariant_vcf,
        popmap = "data/popmap.txt"
    output:
        pi = "results/pixy_output/pi.txt",
        fst =  "results/pixy_output/fst.txt",
        dxy =  "results/pixy_output/dxy.txt"
    params:
        window_size = config["pixy"].get("window_size", 10000),
        output_folder = "results/pixy_output/",
        output_prefix = ""
    conda:
        "../envs/pixy.yaml"
    threads: 24
    resources:
        mem_mb = 8000,
        time = "10:00:00"
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats pi fst dxy --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder}
        """