# Genome scan module - calculates FST, pi, and dXY along the genome
# with configurable window sizes for two populations

rule extract_genome_scan_samples:
    input:
        indpopdata = rules.generate_popdata.output.indpopdata
    output:
        samples_file = "results/{project}/genome_scan/{project}.genome_scan_samples.txt"
    params:
        pop_column = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("population_column", "Group"),
        pop1 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("pop1", ""),
        pop2 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("pop2", "")
    log:
        "logs/{project}/extract_genome_scan_samples.log"
    benchmark:
        "benchmarks/{project}/extract_genome_scan_samples.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/extract_genome_scan_samples.py"

rule subset_vcf_genome_scan:
    input:
        vcf = rules.prepare_invariant_vcf_gz.output.vcf,
        vcf_index = rules.prepare_invariant_vcf_gz_index.output.index,
        samples_file = rules.extract_genome_scan_samples.output.samples_file
    output:
        vcf = "results/{project}/genome_scan/{project}.genome_scan_subset.vcf.gz",
        index = "results/{project}/genome_scan/{project}.genome_scan_subset.vcf.gz.csi"
    log:
        "logs/{project}/subset_vcf_genome_scan.log"
    benchmark:
        "benchmarks/{project}/subset_vcf_genome_scan.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        mkdir -p results/{wildcards.project}/genome_scan
        bcftools view -S {input.samples_file} -Oz -o {output.vcf} {input.vcf} &> {log}
        bcftools index {output.vcf} &> {log}
        """

rule calculate_genome_scan_missing_indv:
    input:
        vcf = rules.subset_vcf_genome_scan.output.vcf
    output:
        imiss = "results/{project}/genome_scan/{project}.genome_scan.imiss"
    log:
        "logs/{project}/calculate_genome_scan_missing_indv.log"
    params:
        out_prefix = "results/{project}/genome_scan/{project}.genome_scan"
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

rule calculate_genome_scan_missing_loci:
    input:
        vcf = rules.subset_vcf_genome_scan.output.vcf
    output:
        lmiss = "results/{project}/genome_scan/{project}.genome_scan.lmiss"
    log:
        "logs/{project}/calculate_genome_scan_missing_loci.log"
    params:
        out_prefix = "results/{project}/genome_scan/{project}.genome_scan"
    conda:
        "../envs/vcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --missing-site \
                 --out {params.out_prefix} &> {log}
        """

rule create_genome_scan_popmap:
    input:
        vcf = rules.subset_vcf_genome_scan.output.vcf,
        indpopdata = rules.generate_popdata.output.indpopdata
    output:
        popmap = "results/{project}/genome_scan/{project}.genome_scan_popmap.txt"
    params:
        pop_column = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("population_column", "Group"),
        pop1 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("pop1", ""),
        pop2 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("pop2", "")
    log:
        "logs/{project}/create_genome_scan_popmap.log"
    benchmark:
        "benchmarks/{project}/create_genome_scan_popmap.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/create_genome_scan_popmap.py"

rule genome_scan_fst:
    input:
        vcf = rules.subset_vcf_genome_scan.output.vcf,
        popmap = rules.create_genome_scan_popmap.output.popmap
    output:
        fst = "results/{project}/genome_scan/{project}.genome_scan_fst.txt"
    log:
        "logs/{project}/genome_scan_fst.log"
    benchmark:
        "benchmarks/{project}/genome_scan_fst.txt"
    params:
        window_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("window_size_fst", 100),
        output_folder = "results/{project}/genome_scan/",
        output_prefix = "{project}.genome_scan"
    conda:
        "../envs/pixy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["genome_scan"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["genome_scan"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["genome_scan"]["runtime"]
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats fst --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder} --output_prefix {params.output_prefix} &> {log}
        """

rule genome_scan_pi:
    input:
        vcf = rules.subset_vcf_genome_scan.output.vcf,
        popmap = rules.create_genome_scan_popmap.output.popmap
    output:
        pi = "results/{project}/genome_scan/{project}.genome_scan_pi.txt"
    log:
        "logs/{project}/genome_scan_pi.log"
    benchmark:
        "benchmarks/{project}/genome_scan_pi.txt"
    params:
        window_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("window_size_pi", 1000000),
        output_folder = "results/{project}/genome_scan/",
        output_prefix = "{project}.genome_scan"
    conda:
        "../envs/pixy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["genome_scan"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["genome_scan"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["genome_scan"]["runtime"]
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats pi --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder} --output_prefix {params.output_prefix} &> {log}
        """

rule genome_scan_dxy:
    input:
        vcf = rules.subset_vcf_genome_scan.output.vcf,
        popmap = rules.create_genome_scan_popmap.output.popmap
    output:
        dxy = "results/{project}/genome_scan/{project}.genome_scan_dxy.txt"
    log:
        "logs/{project}/genome_scan_dxy.log"
    benchmark:
        "benchmarks/{project}/genome_scan_dxy.txt"
    params:
        window_size = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("window_size_dxy", 1000000),
        output_folder = "results/{project}/genome_scan/",
        output_prefix = "{project}.genome_scan"
    conda:
        "../envs/pixy.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["genome_scan"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["genome_scan"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["genome_scan"]["runtime"]
    shell:
        """
        mkdir -p {params.output_folder}
        pixy --stats dxy --vcf {input.vcf} --populations {input.popmap} \
             --n_cores {threads} --window_size {params.window_size} \
             --output_folder {params.output_folder} --output_prefix {params.output_prefix} &> {log}
        """

# Plot missing data histograms
rule plot_genome_scan_imiss_histogram:
    input:
        imiss = rules.calculate_genome_scan_missing_indv.output.imiss
    output:
        pdf = "results/{project}/genome_scan/plots/{project}.genome_scan.imiss_histogram.pdf",
        rds = "results/{project}/genome_scan/plots/{project}.genome_scan.imiss_histogram.rds",
        summary = "results/{project}/genome_scan/{project}.genome_scan.imiss_summary.txt"
    log:
        "logs/{project}/plot_genome_scan_imiss_histogram.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

rule plot_genome_scan_lmiss_histogram:
    input:
        lmiss = rules.calculate_genome_scan_missing_loci.output.lmiss
    output:
        pdf = "results/{project}/genome_scan/plots/{project}.genome_scan.lmiss_histogram.pdf",
        rds = "results/{project}/genome_scan/plots/{project}.genome_scan.lmiss_histogram.rds",
        summary = "results/{project}/genome_scan/{project}.genome_scan.lmiss_summary.txt"
    log:
        "logs/{project}/plot_genome_scan_lmiss_histogram.log"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_missingness_histogram.R"

rule plot_genome_scan:
    input:
        fst = rules.genome_scan_fst.output.fst,
        pi = rules.genome_scan_pi.output.pi,
        dxy = rules.genome_scan_dxy.output.dxy
    output:
        pdf = "results/{project}/genome_scan/plots/{project}.genome_scan_genome_plot.pdf",
        rds = "results/{project}/genome_scan/plots/{project}.genome_scan_genome_plot.rds"
    params:
        pop1 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("pop1", ""),
        pop2 = lambda wildcards: config["projects"][wildcards.project]["parameters"]["genome_scan"].get("pop2", "")
    log:
        "logs/{project}/plot_genome_scan.log"
    benchmark:
        "benchmarks/{project}/plot_genome_scan.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_genome_scan.R"

