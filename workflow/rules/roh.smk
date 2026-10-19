# Rule to index subset VCF
rule bcftools_roh:
    input:
        vcf=rules.subset_vcf_after_relatedness.output.vcf
    output:
        roh="results/{project}/roh/{project}.bcftools_roh.txt"
    log:
        "logs/{project}/bcftools_roh.log"
    benchmark:
        "benchmarks/{project}/bcftools_roh.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    shell:
        """
        set -euo pipefail
        bcftools +fill-tags {input.vcf} -- -t AF | \
        bcftools roh --AF-dflt 0.5 -G 30 -Or -o {output.roh} &> {log}
        """

# 1 Mb windows on the ROH VCF, empty windows included, plus SNPs/Mb inside called ROHs.
rule roh_snp_density:
    input:
        vcf=rules.subset_vcf_after_relatedness.output.vcf,
        roh=rules.bcftools_roh.output.roh
    output:
        windows="results/{project}/roh/{project}.snp_density_windows.txt",
        summary="results/{project}/roh/{project}.snp_density_summary.txt"
    log:
        "logs/{project}/roh_snp_density.log"
    benchmark:
        "benchmarks/{project}/roh_snp_density.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    shell:
        """
        set -euo pipefail
        exec 2> {log}
        {{
          bcftools view -h {input.vcf} | awk '/^##contig=/ {{
            id = ""; len = "NA"
            if (match($0, /ID=[^,>]+/)) id = substr($0, RSTART+3, RLENGTH-3)
            if (match($0, /length=[0-9]+/)) len = substr($0, RSTART+7, RLENGTH-7)
            if (id != "") print "C\\t" id "\\t" len
          }}'
          bcftools query -f '%CHROM\\t%POS\\n' {input.vcf} | awk 'BEGIN {{OFS="\\t"}} {{print "S", $1, $2}}'
          awk -F '\\t' '$1 == "RG" {{print "R\\t" $6 "\\t" $7}}' {input.roh}
        }} | awk -v WIN=1000000 -v windows="{output.windows}" -v summary="{output.summary}" '
        BEGIN {{
            FS = "\\t"
            OFS = "\\t"
            if (WIN == "") WIN = 1000000
        }}
        $1 == "C" {{
            chrom = $2
            if (!(chrom in seen)) {{
                seen[chrom] = 1
                order[++norder] = chrom
            }}
            if ($3 != "NA" && $3 + 0 > 0) length_bp[chrom] = $3 + 0
            next
        }}
        $1 == "S" {{
            chrom = $2
            pos = $3 + 0
            if (pos < 1) next
            if (!(chrom in seen)) {{
                seen[chrom] = 1
                order[++norder] = chrom
            }}
            bin = int((pos - 1) / WIN)
            wcount[chrom SUBSEP bin]++
            nsnp[chrom]++
            if (pos > maxpos[chrom]) maxpos[chrom] = pos
            next
        }}
        $1 == "R" {{
            roh_bp += $2 + 0
            roh_markers += $3 + 0
            next
        }}
        function median_of(a, n,    i, copy) {{
            if (n <= 0) return "NA"
            for (i = 1; i <= n; i++) copy[i] = a[i]
            asort(copy)
            if (n % 2) return copy[(n + 1) / 2]
            return (copy[n / 2] + copy[n / 2 + 1]) / 2
        }}
        END {{
            nchrom = 0
            for (ci = 1; ci <= norder; ci++) {{
                chrom = order[ci]
                span = 0
                if ((chrom in length_bp) && length_bp[chrom] > 0) span = length_bp[chrom]
                if (maxpos[chrom] > span) span = maxpos[chrom]
                if (span < 1) continue
                nchrom++
                chroms[nchrom] = chrom
                spans[nchrom] = span
            }}
            comparable = "no"
            if (nchrom > 0) {{
                med_span = median_of(spans, nchrom)
                if (med_span + 0 >= WIN) comparable = "yes"
            }}

            print "chrom", "window_start", "n_snps" > windows

            g_snps = 0
            g_span = 0
            g_windows = 0
            g_empty = 0
            nall = 0

            for (ci = 1; ci <= nchrom; ci++) {{
                chrom = chroms[ci]
                span = spans[ci]
                nwin = int((span - 1) / WIN) + 1
                delete wc
                empty = 0
                occupied = 0
                snps = nsnp[chrom] + 0
                for (b = 0; b < nwin; b++) {{
                    c = wcount[chrom SUBSEP b] + 0
                    wc[b + 1] = c
                    allc[++nall] = c
                    if (c == 0) empty++
                    else occupied++
                    print chrom, b * WIN, c > windows
                }}
                c_mean = (nwin > 0) ? snps / nwin : "NA"
                c_median = median_of(wc, nwin)
                c_frac = (nwin > 0) ? empty / nwin : "NA"
                c_per_mb = (span > 0) ? snps / (span / 1e6) : "NA"
                c_occ = (occupied > 0) ? snps / occupied : "NA"
                row_chrom[ci] = chrom
                row_snps[ci] = snps
                row_span[ci] = span
                row_nwin[ci] = nwin
                row_mean[ci] = c_mean
                row_median[ci] = c_median
                row_frac[ci] = c_frac
                row_per_mb[ci] = c_per_mb
                row_occ[ci] = c_occ
                g_snps += snps
                g_span += span
                g_windows += nwin
                g_empty += empty
            }}
            close(windows)

            g_mean = (g_windows > 0) ? g_snps / g_windows : "NA"
            g_median = median_of(allc, nall)
            g_frac = (g_windows > 0) ? g_empty / g_windows : "NA"
            g_per_mb = (g_span > 0) ? g_snps / (g_span / 1e6) : "NA"
            g_occupied = g_windows - g_empty
            g_occ = (g_occupied > 0) ? g_snps / g_occupied : "NA"

            print "scope", "chrom", "n_snps", "span_bp", "n_windows", "mean_snps_per_window", "median_snps_per_window", "frac_windows_empty", "snps_per_mb_of_span", "mean_snps_per_occupied_window", "windows_comparable" > summary
            printf "genome\\t.\\t%d\\t%d\\t%d\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", g_snps, g_span, g_windows, num(g_mean), num(g_median), num(g_frac), num(g_per_mb), num(g_occ), comparable > summary
            for (ci = 1; ci <= nchrom; ci++) {{
                printf "chrom\\t%s\\t%d\\t%d\\t%d\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", row_chrom[ci], row_snps[ci], row_span[ci], row_nwin[ci], num(row_mean[ci]), num(row_median[ci]), num(row_frac[ci]), num(row_per_mb[ci]), num(row_occ[ci]), comparable > summary
            }}
            if (roh_bp > 0) {{
                roh_per_mb = roh_markers / (roh_bp / 1e6)
                printf "inside_roh\\t.\\t%d\\t%d\\tNA\\tNA\\tNA\\tNA\\t%s\\tNA\\t%s\\n", roh_markers, roh_bp, num(roh_per_mb), comparable > summary
            }} else {{
                printf "inside_roh\\t.\\t0\\t0\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\t%s\\n", comparable > summary
            }}
            close(summary)
        }}
        function num(x) {{
            if (x == "NA") return "NA"
            return sprintf("%.6f", x + 0)
        }}
        '
        """

rule roh_summary:
    input:
        roh=rules.bcftools_roh.output.roh,
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        summary="results/{project}/roh/{project}.roh_summary.txt",
        per_ind="results/{project}/roh/{project}.roh_per_ind.txt",
        stats="results/{project}/roh/{project}.roh_stats_comparisons.txt"
    params:
        group_by = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("roh", {}).get("group_by", ["Site"])
    log:
        "logs/{project}/{project}.roh_summary.log"
    benchmark:
        "benchmarks/{project}/roh_summary.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/roh_summary.R"


rule roh_plots:
    input:
        roh=rules.bcftools_roh.output.roh,
        per_ind=rules.roh_summary.output.per_ind
    output:
        froh_histogram="results/{project}/roh/plots/{project}.froh_histogram.pdf",
        froh_histogram_rds="results/{project}/roh/plots/{project}.froh_histogram.rds",
        nroh_histogram="results/{project}/roh/plots/{project}.nroh_histogram.pdf",
        nroh_histogram_rds="results/{project}/roh/plots/{project}.nroh_histogram.rds",
        length_histogram="results/{project}/roh/plots/{project}.length_histogram.pdf",
        length_histogram_rds="results/{project}/roh/plots/{project}.length_histogram.rds",
        length_segments_histogram="results/{project}/roh/plots/{project}.length_segments_histogram.pdf",
        length_segments_histogram_rds="results/{project}/roh/plots/{project}.length_segments_histogram.rds",
        froh_classes="results/{project}/roh/plots/{project}.froh_classes.pdf",
        froh_classes_rds="results/{project}/roh/plots/{project}.froh_classes.rds",
        froh_classes_histogram="results/{project}/roh/plots/{project}.froh_classes_histogram.pdf",
        froh_classes_histogram_rds="results/{project}/roh/plots/{project}.froh_classes_histogram.rds",
        nroh_classes="results/{project}/roh/plots/{project}.nroh_classes.pdf",
        nroh_classes_rds="results/{project}/roh/plots/{project}.nroh_classes.rds",
        nroh_classes_histogram="results/{project}/roh/plots/{project}.nroh_classes_histogram.pdf",
        nroh_classes_histogram_rds="results/{project}/roh/plots/{project}.nroh_classes_histogram.rds",
        froh_vs_nroh="results/{project}/roh/plots/{project}.froh_vs_nroh.pdf",
        froh_vs_nroh_rds="results/{project}/roh/plots/{project}.froh_vs_nroh.rds"
    params:
        group_col=None,
        width=lambda wildcards: _roh_plot_style_params(wildcards)["width"],
        height=lambda wildcards: _roh_plot_style_params(wildcards)["height"],
        axis_title_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_text_size"],
        point_size=lambda wildcards: _roh_plot_style_params(wildcards)["point_size"],
    log:
        "logs/{project}/{project}.roh_plots.log"
    benchmark:
        "benchmarks/{project}/roh_plots.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/roh_plots.R"


rule roh_group_plots:
    input:
        per_ind=rules.roh_summary.output.per_ind
    output:
        froh="results/{project}/roh/plots/{project}.froh_by_{group_col}.pdf",
        froh_rds="results/{project}/roh/plots/{project}.froh_by_{group_col}.rds",
        froh_classes="results/{project}/roh/plots/{project}.froh_by_{group_col}_classes.pdf",
        froh_classes_rds="results/{project}/roh/plots/{project}.froh_by_{group_col}_classes.rds",
        nroh="results/{project}/roh/plots/{project}.nroh_by_{group_col}.pdf",
        nroh_rds="results/{project}/roh/plots/{project}.nroh_by_{group_col}.rds",
        nroh_classes="results/{project}/roh/plots/{project}.nroh_by_{group_col}_classes.pdf",
        nroh_classes_rds="results/{project}/roh/plots/{project}.nroh_by_{group_col}_classes.rds"
    params:
        group_col=lambda wildcards: wildcards.group_col,
        group_colors=lambda wildcards: _roh_group_setting(
            wildcards.project, wildcards.group_col, "colors"
        ),
        group_sort_by=lambda wildcards: _roh_group_setting(
            wildcards.project, wildcards.group_col, "sort_by"
        ),
        width=lambda wildcards: _roh_plot_style_params(wildcards)["width"],
        height=lambda wildcards: _roh_plot_style_params(wildcards)["height"],
        axis_title_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_text_size"],
        point_size=lambda wildcards: _roh_plot_style_params(wildcards)["point_size"],
    log:
        "logs/{project}/roh_group_plots.{group_col}.log"
    benchmark:
        "benchmarks/{project}/roh_group_plots.{group_col}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/roh_plots.R"


rule roh_comparison_plots:
    """
    Mean ± SE and Tukey pairwise differences from the ROH comparison report.
    """
    input:
        stats=rules.roh_summary.output.stats
    output:
        pdf="results/{project}/roh/plots/{project}.roh_comparisons_by_{group_col}.pdf",
        rds="results/{project}/roh/plots/{project}.roh_comparisons_by_{group_col}.rds"
    params:
        group_col=lambda wildcards: wildcards.group_col,
        width=lambda wildcards: _roh_plot_style_params(wildcards)["width"],
        height=lambda wildcards: _roh_plot_style_params(wildcards)["height"],
        axis_title_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_title_size"],
        axis_text_size=lambda wildcards: _roh_plot_style_params(wildcards)["axis_text_size"],
        dpi=lambda wildcards: _roh_setting(wildcards.project, "dpi", 300)
    log:
        "logs/{project}/roh_comparison_plots.{group_col}.log"
    benchmark:
        "benchmarks/{project}/roh_comparison_plots.{group_col}.txt"
    conda:
        "../envs/r-plot.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/plot_roh_comparisons.R"


rule plot_froh_map:
    input:
        unpack(_froh_map_inputs),
    output:
        pdf = "results/{project}/roh/plots/{project}.froh_map.pdf",
        rds = "results/{project}/roh/plots/{project}.froh_map.rds"
    params:
        lambda wildcards: _froh_map_params(wildcards)
    log:
        "logs/{project}/plot_froh_map.log"
    benchmark:
        "benchmarks/{project}/plot_froh_map.txt"
    conda:
        "../envs/mapmixture.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["roh"]["runtime"]
    script:
        "../scripts/plot_froh_map.R"
