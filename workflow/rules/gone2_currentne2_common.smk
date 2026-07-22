"""
Shared prep for GONE2 and currentNe2: per-population samples + filtered VCFs.
"""


rule gone2_currentne2_common_prepare_samples:
    input:
        indpopdata=rules.generate_popdata.output.indpopdata,
    output:
        samples_dir=directory("results/{project}/gone2_currentne2_common/samples"),
        populations="results/{project}/gone2_currentne2_common/{project}.gone2_currentne2_common_populations.tsv",
    params:
        population_column=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("population_column", "Site"),
        min_individuals=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("min_individuals", 10),
    log:
        "logs/{project}/gone2_currentne2_common_prepare_samples.log"
    benchmark:
        "benchmarks/{project}/gone2_currentne2_common_prepare_samples.txt"
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"],
    script:
        "../scripts/gone2_currentne2_common_prepare_samples.py"


rule gone2_currentne2_common_subset_vcf:
    input:
        vcf=rules.select_biallelic_snps.output.biallelic_vcf,
        samples=rules.gone2_currentne2_common_prepare_samples.output.samples_dir,
    output:
        vcf="results/{project}/gone2_currentne2_common/vcf/{project}.{stratum}.vcf",
        chrom_filter="results/{project}/gone2_currentne2_common/vcf/{project}.{stratum}.chrom_filter.tsv",
    params:
        samples_file=lambda wildcards: f"results/{wildcards.project}/gone2_currentne2_common/samples/{wildcards.project}.{wildcards.stratum}.samples.txt",
        f_missing=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("f_missing", 1.0),
        mac_threshold=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("mac_threshold", 1),
        min_snps=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("min_snps", 1000),
        recombination_rate=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("recombination_rate_cM_per_Mb", 2.5),
        # GONE2 hard-fails if any chromosome span is <= 20 cM; applied to shared inputs.
        min_chromosome_cM=lambda wildcards: config["projects"][wildcards.project]["parameters"]["gone2_currentne2_common"].get("min_chromosome_cM", 20),
    log:
        "logs/{project}/gone2_currentne2_common_subset_vcf.{stratum}.log"
    benchmark:
        "benchmarks/{project}/gone2_currentne2_common_subset_vcf_{stratum}.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2_currentne2_common", {}).get("threads", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]),
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2_currentne2_common", {}).get("mem_mb", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"]),
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"].get("gone2_currentne2_common", {}).get("runtime", config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.vcf})"
        # -r requires bgzip + index; write compressed temp and index before region subset.
        TMP_VCF="$(mktemp "${{TMPDIR:-/tmp}}/gone2_currentne2_common_subset.XXXXXX.vcf.gz")"
        KEEP_LIST="$(mktemp "${{TMPDIR:-/tmp}}/gone2_currentne2_common_keep.XXXXXX.txt")"
        trap 'rm -f "$TMP_VCF" "$TMP_VCF.tbi" "$TMP_VCF.csi" "$KEEP_LIST"' EXIT

        bcftools view -S {params.samples_file} -Ou {input.vcf} \
        | bcftools filter -e 'F_MISSING > {params.f_missing}' -Ou - \
        | bcftools +fill-tags - -Ou -- -t AC,AN,AF \
        | bcftools view -i 'MAC > {params.mac_threshold}' -Ou - \
        | bcftools view -Oz -o "$TMP_VCF" > {log} 2>&1
        bcftools index --csi "$TMP_VCF" >> {log} 2>&1

        # Drop chromosomes with SNP-span genetic length <= min_chromosome_cM.
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

        echo "LD-Ne chromosome filter written to {output.chrom_filter}" >> {log}
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
