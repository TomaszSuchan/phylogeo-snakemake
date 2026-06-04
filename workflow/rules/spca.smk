# sPCA (Spatial PCA) — adegenet ordination with Moran's eigenvector maps


def _spca_params(wildcards):
    return config["projects"][wildcards.project]["parameters"].get("spca", {})


def _spca_n_axes(wildcards):
    sp = _spca_params(wildcards)
    nfposi = int(sp.get("nfposi", 2))
    nfnega = int(sp.get("nfnega", 2))
    return nfposi + nfnega


rule spca_analysis:
    """
    Run Spatial PCA on the final analysis VCF and sample coordinates from indpopdata.
    Follows adegenet tutorial-spca.pdf (spca, global/local tests, summary).
    """
    input:
        vcf = lambda wildcards: get_filtered_vcf_output(wildcards),
        indpopdata = rules.generate_popdata.output.indpopdata
    output:
        results_rds = "results/{project}/spca/{project}.spca.results.rds",
        scores_global = "results/{project}/spca/{project}.spca.scores.global.txt",
        scores_local = "results/{project}/spca/{project}.spca.scores.local.txt",
        eigenvalues = "results/{project}/spca/{project}.spca.eigenvalues.txt",
        summary_txt = "results/{project}/spca/{project}.spca.summary.txt",
        global_rtest_txt = "results/{project}/spca/{project}.spca.global.rtest.txt",
        local_rtest_txt = "results/{project}/spca/{project}.spca.local.rtest.txt",
        rtest_rds = "results/{project}/spca/{project}.spca.rtest.rds"
    log:
        "logs/{project}/spca_analysis.log"
    benchmark:
        "benchmarks/{project}/spca_analysis.txt"
    params:
        n_pca = lambda wildcards: _spca_params(wildcards).get("n_pca", None),
        nfposi = lambda wildcards: _spca_params(wildcards).get("nfposi", 2),
        nfnega = lambda wildcards: _spca_params(wildcards).get("nfnega", 2),
        type = lambda wildcards: _spca_params(wildcards).get("type", 1),
        d1 = lambda wildcards: _spca_params(wildcards).get("d1", None),
        d2 = lambda wildcards: _spca_params(wildcards).get("d2", None),
        k = lambda wildcards: _spca_params(wildcards).get("k", None),
        a = lambda wildcards: _spca_params(wildcards).get("a", None),
        dmin = lambda wildcards: _spca_params(wildcards).get("dmin", None),
        nperm = lambda wildcards: _spca_params(wildcards).get("nperm", 999),
        rtest_k = lambda wildcards: _spca_params(wildcards).get("rtest_k", 1)
    conda:
        "../envs/spca.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["spca"]["threads"]
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["spca"]["mem_mb"],
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["spca"]["runtime"]
    script:
        "../scripts/spca_analysis.R"


rule spca_plot:
    """
    Core sPCA diagnostic plots (tutorial-spca.pdf):
    barplot(eig), screeplot, plot/plot.spca, colorplot, global/local test histograms, scatters.
    """
    input:
        results_rds = rules.spca_analysis.output.results_rds,
        rtest_rds = rules.spca_analysis.output.rtest_rds,
        indpopdata = rules.generate_popdata.output.indpopdata
    output:
        eig_barplot = "results/{project}/spca/plots/{project}.spca.eig.barplot.pdf",
        screeplot = "results/{project}/spca/plots/{project}.spca.screeplot.pdf",
        composite_plot = "results/{project}/spca/plots/{project}.spca.composite.pdf",
        composite_plot_rds = "results/{project}/spca/plots/{project}.spca.composite.rds",
        map_plot = "results/{project}/spca/plots/{project}.spca.map.pdf",
        map_plot_rds = "results/{project}/spca/plots/{project}.spca.map.rds",
        colorplot_global = "results/{project}/spca/plots/{project}.spca.colorplot.global.pdf",
        colorplot_local = "results/{project}/spca/plots/{project}.spca.colorplot.local.pdf",
        global_rtest_plot = "results/{project}/spca/plots/{project}.spca.global.rtest.pdf",
        local_rtest_plot = "results/{project}/spca/plots/{project}.spca.local.rtest.pdf",
        scatter_global = "results/{project}/spca/plots/{project}.spca.scatter.global.pdf",
        scatter_global_rds = "results/{project}/spca/plots/{project}.spca.scatter.global.rds",
        scatter_local = "results/{project}/spca/plots/{project}.spca.scatter.local.pdf",
        scatter_local_rds = "results/{project}/spca/plots/{project}.spca.scatter.local.rds"
    log:
        "logs/{project}/spca_plot.log"
    benchmark:
        "benchmarks/{project}/spca_plot.txt"
    params:
        type = lambda wildcards: _spca_params(wildcards).get("type", 1),
        color_by = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("spca_plot", {}).get("color_by", "Site"),
        pc_max = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("spca_plot", {}).get("pc_max", 4),
        plot_axis = 0,
        plot_loading = 0
    conda:
        "../envs/spca.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"].get("mem_mb", 8000),
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"].get("runtime", 30)
    script:
        "../scripts/spca_plot.R"


rule spca_axis_plot:
    """
    Per-axis plot(spca, axis=N) panels (tutorial §1.4).
    """
    input:
        results_rds = rules.spca_analysis.output.results_rds
    output:
        axis_plot = "results/{project}/spca/plots/{project}.spca.axis{axis}.pdf"
    log:
        "logs/{project}/spca_axis_plot.{axis}.log"
    benchmark:
        "benchmarks/{project}/spca_axis_plot.{axis}.txt"
    params:
        plot_axis = lambda wildcards: int(wildcards.axis),
        plot_loading = 0,
        type = lambda wildcards: _spca_params(wildcards).get("type", 1),
        color_by = "none",
        pc_max = 2,
        loading_threshold = None
    conda:
        "../envs/spca.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"].get("mem_mb", 4000),
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"].get("runtime", 10)
    script:
        "../scripts/spca_plot.R"


rule spca_loading_plot:
    """
    Per-axis loadingplot($c1^2) allele contributions (tutorial §1.4).
    """
    input:
        results_rds = rules.spca_analysis.output.results_rds
    output:
        loading_plot = "results/{project}/spca/plots/{project}.spca.loading.axis{axis}.pdf"
    log:
        "logs/{project}/spca_loading_plot.{axis}.log"
    benchmark:
        "benchmarks/{project}/spca_loading_plot.{axis}.txt"
    params:
        plot_axis = 0,
        plot_loading = lambda wildcards: int(wildcards.axis),
        type = lambda wildcards: _spca_params(wildcards).get("type", 1),
        color_by = "none",
        pc_max = 2,
        loading_threshold = lambda wildcards: config["projects"][wildcards.project]["parameters"].get("spca_plot", {}).get("loading_threshold", None)
    conda:
        "../envs/spca.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"].get("mem_mb", 4000),
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"].get("runtime", 10)
    script:
        "../scripts/spca_plot.R"
