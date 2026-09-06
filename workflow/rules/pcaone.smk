# Rule for PCAone analysis using EMU
rule pcaone_emu:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = "results/{project}/pcaone_EMU/{project}.PCA_EMU.eigvecs",
        pcaone_eigenvectors2 = "results/{project}/pcaone_EMU/{project}.PCA_EMU.eigvecs2",
        pcaone_eigenvalues = "results/{project}/pcaone_EMU/{project}.PCA_EMU.eigvals"
    log:
        "logs/{project}/pcaone_emu.log"
    benchmark:
        "benchmarks/{project}/pcaone_emu.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("EMU_SVD_method", 2),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = "results/{project}/pcaone_EMU/{project}.PCA_EMU",
        # Get the bfile prefix (remove .bed extension)
        bfile_prefix = lambda wildcards, input: input.bed.replace('.bed', '')
    conda:
        "../envs/pcaone.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone_emu"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone_emu"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone_emu"]["runtime"]
    shell:
        """
        PCAone --threads {threads} \
        -d {params.SVD_method} \
        --pc {params.PCnum} \
        --emu \
        --bfile {params.bfile_prefix} \
        --out {params.output_prefix} &> {log}
        """

# Rule to run PCAone analysis
rule pcaone:
    input:
        bed = rules.vcf_to_plink.output.bed,
        bim = rules.vcf_to_plink.output.bim,
        fam = rules.vcf_to_plink.output.fam
    output:
        pcaone_eigenvectors = "results/{project}/pcaone/{project}.PCA.eigvecs",
        pcaone_eigenvectors2 = "results/{project}/pcaone/{project}.PCA.eigvecs2",
        pcaone_eigenvalues = "results/{project}/pcaone/{project}.PCA.eigvals",
        pcaone_loadings = "results/{project}/pcaone/{project}.PCA.loadings",
        pcaone_mbim = "results/{project}/pcaone/{project}.PCA.mbim"
    log:
        "logs/{project}/pcaone.log"
    benchmark:
        "benchmarks/{project}/pcaone.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("SVD_method", 3),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = "results/{project}/pcaone/{project}.PCA",
        # Get the bfile prefix (remove .bed extension)
        bfile_prefix = lambda wildcards, input: input.bed.replace('.bed', '')
    conda:
        "../envs/pcaone.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["runtime"]
    shell:
        """
        PCAone --threads {threads} \
        -d {params.SVD_method} \
        --pc {params.PCnum} \
        --printv \
        --bfile {params.bfile_prefix} \
        --out {params.output_prefix} &> {log}
        """

# Rule to build a per-SNP PCA loadings table (which SNPs drive each axis).
# PCAone's .loadings file has one row per SNP (same order as the PLINK .bed,
# which preserves the input VCF order) and one column per PC. We attach
# snp_id/chrom/pos from the *VCF* rather than the .bim, because the PLINK export
# uses `--allow-extra-chr 0` and flattens every non-integer contig name to "0".
# The VCF keeps the real chromosome (e.g. LR999934.1) and position.
rule pcaone_loadings_table:
    input:
        loadings = rules.pcaone.output.pcaone_loadings,
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards)
    output:
        table = "results/{project}/pcaone/{project}.PCA.loadings.tsv"
    log:
        "logs/{project}/pcaone_loadings_table.log"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        r"""
        (
        nvcf=$(zcat {input.vcf} | grep -vc '^#')
        nload=$(wc -l < {input.loadings})
        if [ "$nvcf" -ne "$nload" ]; then
            echo "ERROR: VCF variant count ($nvcf) and .loadings ($nload) row counts differ; cannot join by order" >&2
            exit 1
        fi
        paste \
            <(zcat {input.vcf} | awk -F'\t' 'BEGIN{{OFS="\t"}} !/^#/ {{print $3, $1, $2}}') \
            <(awk 'BEGIN{{OFS="\t"}} {{$1=$1; print}}' {input.loadings}) \
        | awk 'BEGIN{{OFS="\t"}}
               NR==1 {{
                   printf "snp_id\tchrom\tpos"
                   for (i=4; i<=NF; i++) printf "\tPC%d", i-3
                   printf "\n"
               }}
               {{print}}'
        ) > {output.table} 2> {log}
        """


# Rule to summarise which SNPs drive each PCA axis. Produces two tables from the
# per-SNP loadings, restricted to the leading `n_pc` PCs (the rest are noise):
#   * every SNP ranked by |loading| within each PC  (long/tidy)
#   * every SNP assigned to its dominant axis (PC with the largest |loading|)
# Both carry the real chrom/pos recovered in pcaone_loadings_table.
rule pcaone_top_loadings:
    input:
        table = rules.pcaone_loadings_table.output.table
    output:
        ranked = "results/{project}/pcaone/{project}.PCA.loadings.ranked_per_PC.tsv",
        assign = "results/{project}/pcaone/{project}.PCA.loadings.axis_assignment.tsv"
    params:
        # Summarise the same leading PCs that were computed (PCnum); the trailing
        # PCs are noise and are intentionally not summarised.
        n_pc = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10)
    log:
        "logs/{project}/pcaone_top_loadings.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        r"""
        (
        # Every SNP ranked by absolute loading within each PC. loading_sq
        # (= loading^2) is the fraction of that PC's variance carried by the SNP
        # (they sum to 1 across all SNPs, so the per-SNP average is 1/M).
        tail -n +2 {input.table} \
        | awk -F'\t' -v K={params.n_pc} 'BEGIN{{OFS="\t"}}
               {{ for (k=1; k<=K; k++) {{ c=3+k; v=$c; a=(v<0)?-v:v; print k, a, $1, $2, $3, v }} }}' \
        | LC_ALL=C sort -k1,1n -k2,2gr \
        | awk -F'\t' 'BEGIN{{OFS="\t"; print "pc","rank","snp_id","chrom","pos","loading","abs_loading","loading_sq"}}
               {{ cnt[$1]++; print "PC"$1, cnt[$1], $3, $4, $5, $6, $2, $2*$2 }}' \
        > {output.ranked}

        # Assign every SNP to the PC (among the leading K) where |loading| is largest.
        tail -n +2 {input.table} \
        | awk -F'\t' -v K={params.n_pc} 'BEGIN{{OFS="\t"; print "snp_id","chrom","pos","best_pc","loading","abs_loading","loading_sq"}}
               {{ best=0; bestv=-1; bestl=0;
                  for (k=1; k<=K; k++) {{ c=3+k; v=$c; a=(v<0)?-v:v; if (a>bestv) {{ bestv=a; best=k; bestl=v }} }}
                  print $1, $2, $3, "PC"best, bestl, bestv, bestv*bestv }}' \
        > {output.assign}
        ) 2> {log}
        """

# Rule to run PCAone for each miss data threshold
rule pcaone_miss:
    input:
        bed = rules.missing_vcf_to_plink.output.bed,
        bim = rules.missing_vcf_to_plink.output.bim,
        fam = rules.missing_vcf_to_plink.output.fam
    output:
        eigenvectors = "results/{project}/pcaone_miss{miss}/{project}.PCA_miss{miss}.eigvecs",
        eigenvectors2 = "results/{project}/pcaone_miss{miss}/{project}.PCA_miss{miss}.eigvecs2",
        eigenvalues = "results/{project}/pcaone_miss{miss}/{project}.PCA_miss{miss}.eigvals"
    log:
        "logs/{project}/pcaone_miss_{miss}.log"
    benchmark:
        "benchmarks/{project}/pcaone_miss_{miss}.txt"
    params:
        SVD_method = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("SVD_method", 3),
        PCnum = lambda wildcards: config["projects"][wildcards.project]["parameters"]["PCAone"].get("PCnum", 10),
        output_prefix = lambda wildcards: f"results/{wildcards.project}/pcaone_miss{wildcards.miss}/{wildcards.project}.PCA_miss{wildcards.miss}",
        bfile_prefix = lambda wildcards: f"results/{wildcards.project}/filtered_data/{wildcards.project}.biallelic_snps_thinned_miss{wildcards.miss}"
    conda:
        "../envs/pcaone.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["pcaone"]["runtime"]
    shell:
        """
        PCAone --threads {threads} \
        -d {params.SVD_method} \
        --pc {params.PCnum} \
        --bfile {params.bfile_prefix} \
        --out {params.output_prefix} &> {log}
        """

# Rule to plot PCA (colored by population/metadata)
rule plot_pca_colored:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-{color_by}.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_colored_PC{pc1}_PC{pc2}_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        color_by = lambda wildcards: wildcards.color_by,
        group_colors = lambda wildcards: _pca_plot_group_setting(
            wildcards.project, wildcards.color_by, "colors"
        ),
        plot_type = "colored",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA with labels only
rule plot_pca_labeled:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-labeled.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-labeled.rds"
    log:
        "logs/{project}/plot_pca_labeled_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "labeled",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA colored by missing data
rule plot_pca_missing:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-missing.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-PC{pc1}_PC{pc2}-missing.rds"
    log:
        "logs/{project}/plot_pca_missing_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "missing",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA facet (all PC combinations, colored)
rule plot_pca_facet_colored:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-facet-{color_by}.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-facet-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_facet_colored_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        group_colors = lambda wildcards: _pca_plot_group_setting(
            wildcards.project, wildcards.color_by, "colors"
        ),
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "colored",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA facet (all PC combinations, labeled)
rule plot_pca_facet_labeled:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-facet-labeled.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-facet-labeled.rds"
    log:
        "logs/{project}/plot_pca_facet_labeled.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "labeled",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA facet (all PC combinations, missing)
rule plot_pca_facet_missing:
    input:
        eigvecs=rules.pcaone.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone/plots/{project}.PCA-facet-missing.pdf",
        rds="results/{project}/pcaone/plots/{project}.PCA-facet-missing.rds"
    log:
        "logs/{project}/plot_pca_facet_missing.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "missing",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA EMU (colored by population/metadata)
rule plot_pca_emu_colored:
    input:
        eigvecs=rules.pcaone_emu.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone_emu.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-{color_by}.pdf",
        rds="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_emu_colored_PC{pc1}_PC{pc2}_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        color_by = lambda wildcards: wildcards.color_by,
        group_colors = lambda wildcards: _pca_plot_group_setting(
            wildcards.project, wildcards.color_by, "colors"
        ),
        plot_type = "colored",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA EMU with labels only
rule plot_pca_emu_labeled:
    input:
        eigvecs=rules.pcaone_emu.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone_emu.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-labeled.pdf",
        rds="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-labeled.rds"
    log:
        "logs/{project}/plot_pca_emu_labeled_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "labeled",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA EMU colored by missing data
rule plot_pca_emu_missing:
    input:
        eigvecs=rules.pcaone_emu.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone_emu.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-missing.pdf",
        rds="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-PC{pc1}_PC{pc2}-missing.rds"
    log:
        "logs/{project}/plot_pca_emu_missing_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "missing",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA EMU facet (all PC combinations, colored)
rule plot_pca_emu_facet_colored:
    input:
        eigvecs=rules.pcaone_emu.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone_emu.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-facet-{color_by}.pdf",
        rds="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-facet-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_emu_facet_colored_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        group_colors = lambda wildcards: _pca_plot_group_setting(
            wildcards.project, wildcards.color_by, "colors"
        ),
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "colored",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA EMU facet (all PC combinations, labeled)
rule plot_pca_emu_facet_labeled:
    input:
        eigvecs=rules.pcaone_emu.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone_emu.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-facet-labeled.pdf",
        rds="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-facet-labeled.rds"
    log:
        "logs/{project}/plot_pca_emu_facet_labeled.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "labeled",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA EMU facet (all PC combinations, missing)
rule plot_pca_emu_facet_missing:
    input:
        eigvecs=rules.pcaone_emu.output.pcaone_eigenvectors2,
        eigvals=rules.pcaone_emu.output.pcaone_eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_thinned.output.imiss
    output:
        pdf="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-facet-missing.pdf",
        rds="results/{project}/pcaone_EMU/plots/{project}.PCA_EMU-facet-missing.rds"
    log:
        "logs/{project}/plot_pca_emu_facet_missing.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "missing",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA for each miss threshold (colored by population/metadata)
rule plot_pca_miss_colored:
    input:
        eigvecs=rules.pcaone_miss.output.eigenvectors2,
        eigvals=rules.pcaone_miss.output.eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_miss.output.imiss
    output:
        pdf="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-{color_by}.pdf",
        rds="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_miss{miss}_colored_PC{pc1}_PC{pc2}_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        color_by = lambda wildcards: wildcards.color_by,
        group_colors = lambda wildcards: _pca_plot_group_setting(
            wildcards.project, wildcards.color_by, "colors"
        ),
        plot_type = "colored",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA for each miss threshold with labels only
rule plot_pca_miss_labeled:
    input:
        eigvecs=rules.pcaone_miss.output.eigenvectors2,
        eigvals=rules.pcaone_miss.output.eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_miss.output.imiss
    output:
        pdf="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-labeled.pdf",
        rds="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-labeled.rds"
    log:
        "logs/{project}/plot_pca_miss{miss}_labeled_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "labeled",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA for each miss threshold colored by missing data
rule plot_pca_miss_missing:
    input:
        eigvecs=rules.pcaone_miss.output.eigenvectors2,
        eigvals=rules.pcaone_miss.output.eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_miss.output.imiss
    output:
        pdf="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-missing.pdf",
        rds="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-PC{pc1}_PC{pc2}-missing.rds"
    log:
        "logs/{project}/plot_pca_miss{miss}_missing_PC{pc1}_PC{pc2}.log"
    params:
        pc1 = lambda wildcards: wildcards.pc1,
        pc2 = lambda wildcards: wildcards.pc2,
        plot_type = "missing",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_single.R"

# Rule to plot PCA facet for each miss threshold (all PC combinations, colored)
rule plot_pca_miss_facet_colored:
    input:
        eigvecs=rules.pcaone_miss.output.eigenvectors2,
        eigvals=rules.pcaone_miss.output.eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_miss.output.imiss
    output:
        pdf="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-facet-{color_by}.pdf",
        rds="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-facet-{color_by}.rds"
    log:
        "logs/{project}/plot_pca_miss{miss}_facet_colored_{color_by}.log"
    wildcard_constraints:
        color_by="(?!labeled|missing).*"
    params:
        color_by = lambda wildcards: wildcards.color_by,
        group_colors = lambda wildcards: _pca_plot_group_setting(
            wildcards.project, wildcards.color_by, "colors"
        ),
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "colored",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA facet for each miss threshold (all PC combinations, labeled)
rule plot_pca_miss_facet_labeled:
    input:
        eigvecs=rules.pcaone_miss.output.eigenvectors2,
        eigvals=rules.pcaone_miss.output.eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_miss.output.imiss
    output:
        pdf="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-facet-labeled.pdf",
        rds="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-facet-labeled.rds"
    log:
        "logs/{project}/plot_pca_miss{miss}_facet_labeled.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "labeled",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"

# Rule to plot PCA facet for each miss threshold (all PC combinations, missing)
rule plot_pca_miss_facet_missing:
    input:
        eigvecs=rules.pcaone_miss.output.eigenvectors2,
        eigvals=rules.pcaone_miss.output.eigenvalues,
        indpopdata=rules.generate_popdata.output.indpopdata,
        indmiss=rules.calculate_missing_indv_miss.output.imiss
    output:
        pdf="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-facet-missing.pdf",
        rds="results/{project}/pcaone_miss{miss}/plots/{project}.PCA_miss{miss}-facet-missing.rds"
    log:
        "logs/{project}/plot_pca_miss{miss}_facet_missing.log"
    params:
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("pca_plot", {}).get("pc_max", 2),
        plot_type = "missing",
        axis_title_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_title_size", 10),
        axis_text_size = lambda wildcards: _pca_plot_setting(wildcards.project, "axis_text_size", 8),
        point_size = lambda wildcards: _pca_plot_setting(wildcards.project, "point_size", 3),
        width = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "width"), 15.24),
        height = lambda wildcards: _fig_cm_to_in(_pca_plot_setting(wildcards.project, "height"), 12.7),
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    group: "plot_pca"
    conda:
        "../envs/r-plot.yaml"
    script:
        "../scripts/plot_pca_facet.R"
