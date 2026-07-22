"""
GONE2: historical Ne from linkage disequilibrium (https://github.com/esrud/GONE2).
"""


rule gone2_install:
    output:
        sentinel=touch(".snakemake/gone2_installed"),
    params:
        install_dir=".snakemake/gone2",
    log:
        "logs/gone2_install.log"
    conda:
        "../envs/gone2.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        INSTALL_DIR="{params.install_dir}"
        mkdir -p "$INSTALL_DIR"
        if [ ! -x "$INSTALL_DIR/gone2" ]; then
            rm -rf "$INSTALL_DIR/src"
            git clone --depth 1 https://github.com/esrud/GONE2 "$INSTALL_DIR/src" >> {log} 2>&1
            cd "$INSTALL_DIR/src"
            if [ "$(uname -s)" = "Darwin" ]; then
                make macos >> {log} 2>&1 || make gone >> {log} 2>&1
            else
                make gone >> {log} 2>&1
            fi
            test -x gone2 || {{ echo "ERROR: gone2 binary was not built" >> {log}; exit 1; }}
            cp -f gone2 "$INSTALL_DIR/gone2"
            chmod +x "$INSTALL_DIR/gone2"
        fi
        touch {output.sentinel}
        """


rule gone2_prepare_samples:
    input:
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        samples_dir=directory("results/{project}/gone2/samples"),
        populations="results/{project}/gone2/{project}.gone2_populations.tsv",
    params:
        population_column=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("population_column", "Site"),
        min_individuals=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("min_individuals", 10),
    log:
        "logs/{project}/gone2_prepare_samples.log"
    benchmark:
        "benchmarks/{project}/gone2_prepare_samples.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/gone2_prepare_samples.py"


rule gone2_subset_vcf:
    input:
        vcf=rules.select_biallelic_snps.output.biallelic_vcf,
        samples=rules.gone2_prepare_samples.output.samples_dir,
    output:
        vcf="results/{project}/gone2/vcf/{project}.{stratum}.vcf",
        chrom_filter="results/{project}/gone2/vcf/{project}.{stratum}.gone2_chrom_filter.tsv",
    params:
        samples_file=lambda wildcards: f"results/{wildcards.project}/gone2/samples/{wildcards.project}.{wildcards.stratum}.samples.txt",
        f_missing=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("f_missing", 1.0),
        mac_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("mac_threshold", 1),
        min_snps=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("min_snps", 1000),
        recombination_rate=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("recombination_rate_cM_per_Mb", 2.5),
        # GONE2 hard-fails if any chromosome span is <= 20 cM.
        min_chromosome_cM=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("min_chromosome_cM", 20),
    log:
        "logs/{project}/gone2_subset_vcf.{stratum}.log"
    benchmark:
        "benchmarks/{project}/gone2_subset_vcf_{stratum}.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2", {}).get("threads", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]),
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2", {}).get("mem_mb", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"]),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2", {}).get("runtime", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.vcf})"
        # -r requires bgzip + index; write compressed temp and index before region subset.
        TMP_VCF="$(mktemp "${{TMPDIR:-/tmp}}/gone2_subset.XXXXXX.vcf.gz")"
        KEEP_LIST="$(mktemp "${{TMPDIR:-/tmp}}/gone2_keep.XXXXXX.txt")"
        trap 'rm -f "$TMP_VCF" "$TMP_VCF.tbi" "$TMP_VCF.csi" "$KEEP_LIST"' EXIT

        bcftools view -S {params.samples_file} -Ou {input.vcf} \
        | bcftools filter -e 'F_MISSING > {params.f_missing}' -Ou - \
        | bcftools +fill-tags - -Ou -- -t AC,AN,AF \
        | bcftools view -i 'MAC > {params.mac_threshold}' -Ou - \
        | bcftools view -Oz -o "$TMP_VCF" > {log} 2>&1
        bcftools index --csi "$TMP_VCF" >> {log} 2>&1

        # Drop chromosomes with SNP-span genetic length <= min_chromosome_cM (GONE2 requirement).
        bcftools query -f '%CHROM\t%POS\n' "$TMP_VCF" \
        | awk -v rate={params.recombination_rate} -v mincm={params.min_chromosome_cM} \
            -v report={output.chrom_filter:q} -v keep="$KEEP_LIST" '
            BEGIN {{
                OFS = "\t"
                print "chrom", "n_snps", "min_pos", "max_pos", "span_bp", "span_cM", \
                      "recombination_rate_cM_per_Mb", "min_chromosome_cM", "status" > report
            }}
            {{
                c = $1; p = $2 + 0
                if (!(c in min) || p < min[c]) min[c] = p
                if (!(c in max) || p > max[c]) max[c] = p
                n[c]++
            }}
            END {{
                for (c in n) {{
                    span = max[c] - min[c]
                    if (span < 0) span = 0
                    cm = span / 1e6 * rate
                    status = (cm > mincm) ? "KEEP" : "DROP"
                    print c, n[c], min[c], max[c], span, cm, rate, mincm, status >> report
                    if (status == "KEEP") print c >> keep
                }}
            }}'

        echo "GONE2 chromosome filter written to {output.chrom_filter}" >> {log}
        if [ ! -s "$KEEP_LIST" ]; then
            echo "ERROR: no chromosomes longer than {params.min_chromosome_cM} cM at {params.recombination_rate} cM/Mb" >> {log}
            exit 1
        fi
        REGIONS="$(paste -sd, "$KEEP_LIST")"
        echo "Keeping chromosomes: $REGIONS" >> {log}
        bcftools view -r "$REGIONS" -Ov -o {output.vcf} "$TMP_VCF" >> {log} 2>&1

        NSNPS=$(grep -vc '^#' {output.vcf} || true)
        if [ "$NSNPS" -lt {params.min_snps} ]; then
            echo "ERROR: only $NSNPS SNPs after filtering (min_snps={params.min_snps})" >> {log}
            exit 1
        fi
        """

rule gone2_run:
    input:
        install=rules.gone2_install.output.sentinel,
        vcf=rules.gone2_subset_vcf.output.vcf,
        # Forces re-subset when chrom filtering was added after existing VCFs.
        chrom_filter=rules.gone2_subset_vcf.output.chrom_filter,
    output:
        ne="results/{project}/gone2/vcf/{project}.{stratum}_GONE2_Ne",
    params:
        gone2_bin=".snakemake/gone2/gone2",
        # GONE2 appends _GONE2_Ne to -o; keep stem aligned with output path.
        out_stem=lambda wildcards, output: str(output.ne).removesuffix("_GONE2_Ne"),
        recombination_rate=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("recombination_rate_cM_per_Mb", 2.5),
    log:
        "logs/{project}/gone2_run.{stratum}.log"
    benchmark:
        "benchmarks/{project}/gone2_run_{stratum}.txt"
    conda:
        "../envs/gone2.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2", {}).get("threads", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]),
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2", {}).get("mem_mb", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"]),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2", {}).get("runtime", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]),
    shell:
        """
        {params.gone2_bin} -g 0 -r {params.recombination_rate} -t {threads} \
            -o {params.out_stem} {input.vcf} > {log} 2>&1
        """


rule gone2_plot_ne:
    input:
        ne=rules.gone2_run.output.ne,
    output:
        pdf="results/{project}/gone2/plots/{project}.{stratum}.gone2_ne.pdf",
        rds="results/{project}/gone2/plots/{project}.{stratum}.gone2_ne.rds",
    params:
        width=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get("width", 8),
        height=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get("height", 5),
        dpi=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2"].get("plot", {}).get("dpi", 300),
    log:
        "logs/{project}/gone2_plot_ne.{stratum}.log"
    benchmark:
        "benchmarks/{project}/gone2_plot_ne_{stratum}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/plot_gone2_ne.R"


rule gone2_collect_summary:
    input:
        unpack(_gone2_summary_inputs),
    output:
        summary="results/{project}/gone2/{project}.gone2_summary.tsv",
    params:
        vcf_dir=lambda wildcards: f"results/{wildcards.project}/gone2/vcf",
    log:
        "logs/{project}/gone2_collect_summary.log"
    benchmark:
        "benchmarks/{project}/gone2_collect_summary.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/gone2_collect_summary.py"
