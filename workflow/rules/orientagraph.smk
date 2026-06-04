# OrientAGraph historical population graph inference


def _treemix_params(wildcards):
    return config["projects"][wildcards.project]["parameters"].get("treemix", {})


def _treemix_migration_edges_from_params(params):
    tm = params.get("treemix", {})
    edges = tm.get("migration_edges", [0])
    if isinstance(edges, int):
        edges = [edges]
    return [int(edge) for edge in edges]


def _treemix_migration_edges(wildcards):
    return _treemix_migration_edges_from_params(config["projects"][wildcards.project]["parameters"])


def _treemix_optm_replicates_from_params(params):
    tm = params.get("treemix", {})
    optm = tm.get("optm", {})
    return list(range(1, int(optm.get("replicates", 1)) + 1))


def _treemix_optm_replicates(wildcards):
    return _treemix_optm_replicates_from_params(config["projects"][wildcards.project]["parameters"])


def _treemix_bootstrap_config_from_params(params):
    tm = params.get("treemix", {})
    bootstrap = tm.get("bootstrap", {})
    if isinstance(bootstrap, bool):
        return {"enabled": bootstrap, "migration_edges": [], "replicates": 100}
    return bootstrap or {}


def _treemix_bootstrap_enabled_from_params(params):
    return bool(_treemix_bootstrap_config_from_params(params).get("enabled", False))


def _treemix_bootstrap_migration_edges_from_params(params):
    bootstrap = _treemix_bootstrap_config_from_params(params)
    edges = bootstrap.get("migration_edges", [])
    if isinstance(edges, int):
        edges = [edges]
    edges = [int(edge) for edge in edges]
    if _treemix_bootstrap_enabled_from_params(params) and not edges:
        raise ValueError("treemix.bootstrap.migration_edges must be set when treemix.bootstrap.enabled is true")
    return edges


def _treemix_bootstrap_migration_edges(wildcards):
    return _treemix_bootstrap_migration_edges_from_params(config["projects"][wildcards.project]["parameters"])


def _treemix_bootstrap_replicates_from_params(params):
    bootstrap = _treemix_bootstrap_config_from_params(params)
    return list(range(1, int(bootstrap.get("replicates", 100)) + 1))


def _treemix_bootstrap_replicates(wildcards):
    return _treemix_bootstrap_replicates_from_params(config["projects"][wildcards.project]["parameters"])


def _treemix_main_bootstrap_flag(wildcards):
    # Backward compatibility for a legacy boolean; the new bootstrap block is handled by dedicated rules.
    return "-bootstrap" if _treemix_params(wildcards).get("bootstrap", False) is True else ""


def _orientagraph_value_flag(wildcards, key, flag, default=None):
    value = _treemix_params(wildcards).get(key, default)
    if value in [None, False, "", "null", "NULL"]:
        return ""
    if value is True:
        return flag
    value = str(value)
    if value.lower() in {"all", "every", "each"}:
        return flag
    return f"{flag} {value}"


rule treemix_prepare_input:
    """
    Convert the final analysis VCF to TreeMix-compatible allele-count input.
    Population groups are taken from a configured indpopdata column.
    """
    input:
        vcf=lambda wildcards: get_filtered_vcf_output(wildcards),
        indpopdata=rules.generate_popdata.output.indpopdata
    output:
        frq="results/{project}/treemix/{project}.treemix.frq.gz",
        clust="results/{project}/treemix/{project}.treemix.clust",
        populations="results/{project}/treemix/{project}.treemix.populations.tsv",
        positions="results/{project}/treemix/{project}.treemix.positions.tsv",
        summary="results/{project}/treemix/{project}.treemix.input_summary.tsv"
    params:
        population_column=lambda wildcards: _treemix_params(wildcards).get("population_column", "Site"),
        min_called_populations=lambda wildcards: _treemix_params(wildcards).get("min_called_populations", 2)
    log:
        "logs/{project}/treemix_prepare_input.log"
    benchmark:
        "benchmarks/{project}/treemix_prepare_input.txt"
    conda:
        "../envs/python.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/vcf_to_treemix.py"


rule treemix:
    """
    Fit an OrientAGraph graph with a configured number of migration edges.
    """
    input:
        frq=rules.treemix_prepare_input.output.frq
    output:
        llik="results/{project}/treemix/{project}.treemix.m{m}.llik",
        treeout="results/{project}/treemix/{project}.treemix.m{m}.treeout.gz",
        cov="results/{project}/treemix/{project}.treemix.m{m}.cov.gz",
        covse="results/{project}/treemix/{project}.treemix.m{m}.covse.gz",
        modelcov="results/{project}/treemix/{project}.treemix.m{m}.modelcov.gz",
        vertices="results/{project}/treemix/{project}.treemix.m{m}.vertices.gz",
        edges="results/{project}/treemix/{project}.treemix.m{m}.edges.gz"
    params:
        prefix="results/{project}/treemix/{project}.treemix.m{m}",
        block_size=lambda wildcards: _treemix_params(wildcards).get(
            "k", _treemix_params(wildcards).get("block_size", 1000)
        ),
        root_flag=lambda wildcards: (
            "-root " + str(_treemix_params(wildcards).get("root"))
            if _treemix_params(wildcards).get("root") not in [None, "", "null", "NULL"]
            else ""
        ),
        global_flag=lambda wildcards: "-global" if _treemix_params(wildcards).get("global", False) else "",
        bootstrap_flag=lambda wildcards: _treemix_main_bootstrap_flag(wildcards),
        mlno_flag=lambda wildcards: _orientagraph_value_flag(wildcards, "mlno", "-mlno", True),
        allmigs_flag=lambda wildcards: _orientagraph_value_flag(wildcards, "allmigs", "-allmigs", False),
        seed_flag=lambda wildcards: (
            "-seed " + str(_treemix_params(wildcards).get("seed"))
            if _treemix_params(wildcards).get("seed") not in [None, "", "null", "NULL"]
            else ""
        )
    log:
        "logs/{project}/treemix.m{m}.log"
    benchmark:
        "benchmarks/{project}/treemix.m{m}.txt"
    conda:
        "../envs/orientagraph.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["runtime"]
    shell:
        """
        orientagraph \
            -i {input.frq} \
            -o {params.prefix} \
            -m {wildcards.m} \
            -k {params.block_size} \
            {params.root_flag} \
            {params.global_flag} \
            {params.bootstrap_flag} \
            {params.mlno_flag} \
            {params.allmigs_flag} \
            {params.seed_flag} \
            > {log} 2>&1
        """


rule treemix_optm_run:
    """
    Additional OrientAGraph replicate runs used by OptM migration-edge selection.
    """
    input:
        frq=rules.treemix_prepare_input.output.frq
    output:
        llik="results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.llik",
        treeout="results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.treeout.gz",
        cov="results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.cov.gz",
        covse="results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.covse.gz",
        modelcov="results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.modelcov.gz",
        vertices="results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.vertices.gz",
        edges="results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.edges.gz"
    params:
        prefix="results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}",
        block_size=lambda wildcards: _treemix_params(wildcards).get(
            "k", _treemix_params(wildcards).get("block_size", 1000)
        ),
        root_flag=lambda wildcards: (
            "-root " + str(_treemix_params(wildcards).get("root"))
            if _treemix_params(wildcards).get("root") not in [None, "", "null", "NULL"]
            else ""
        ),
        global_flag=lambda wildcards: "-global" if _treemix_params(wildcards).get("global", False) else "",
        bootstrap_flag=lambda wildcards: _treemix_main_bootstrap_flag(wildcards),
        mlno_flag=lambda wildcards: _orientagraph_value_flag(wildcards, "mlno", "-mlno", True),
        allmigs_flag=lambda wildcards: _orientagraph_value_flag(wildcards, "allmigs", "-allmigs", False),
        seed_flag=lambda wildcards: (
            "-seed " + str(
                int(_treemix_params(wildcards).get("seed")) +
                int(wildcards.rep) * 1000 +
                int(wildcards.m)
            )
            if _treemix_params(wildcards).get("seed") not in [None, "", "null", "NULL"]
            else ""
        )
    log:
        "logs/{project}/treemix_optm.r{rep}.m{m}.log"
    benchmark:
        "benchmarks/{project}/treemix_optm.r{rep}.m{m}.txt"
    conda:
        "../envs/orientagraph.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["runtime"]
    shell:
        """
        orientagraph \
            -i {input.frq} \
            -o {params.prefix} \
            -m {wildcards.m} \
            -k {params.block_size} \
            {params.root_flag} \
            {params.global_flag} \
            {params.bootstrap_flag} \
            {params.mlno_flag} \
            {params.allmigs_flag} \
            {params.seed_flag} \
            > {log} 2>&1
        """


rule treemix_bootstrap_run:
    """
    Run one OrientAGraph / TreeMix bootstrap replicate for a selected final migration-edge count.
    """
    input:
        frq=rules.treemix_prepare_input.output.frq
    output:
        llik="results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.llik",
        treeout="results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.treeout.gz",
        cov="results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.cov.gz",
        covse="results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.covse.gz",
        modelcov="results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.modelcov.gz",
        vertices="results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.vertices.gz",
        edges="results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.edges.gz"
    params:
        prefix="results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}",
        block_size=lambda wildcards: _treemix_params(wildcards).get(
            "k", _treemix_params(wildcards).get("block_size", 1000)
        ),
        root_flag=lambda wildcards: (
            "-root " + str(_treemix_params(wildcards).get("root"))
            if _treemix_params(wildcards).get("root") not in [None, "", "null", "NULL"]
            else ""
        ),
        global_flag=lambda wildcards: "-global" if _treemix_params(wildcards).get("global", False) else "",
        mlno_flag=lambda wildcards: _orientagraph_value_flag(wildcards, "mlno", "-mlno", True),
        allmigs_flag=lambda wildcards: _orientagraph_value_flag(wildcards, "allmigs", "-allmigs", False),
        seed_flag=lambda wildcards: "-seed " + str(
            (
                int(_treemix_params(wildcards).get("seed"))
                if _treemix_params(wildcards).get("seed") not in [None, "", "null", "NULL"]
                else 0
            ) +
            int(wildcards.rep) * 10000 +
            int(wildcards.m)
        )
    log:
        "logs/{project}/treemix_bootstrap.r{rep}.m{m}.log"
    benchmark:
        "benchmarks/{project}/treemix_bootstrap.r{rep}.m{m}.txt"
    conda:
        "../envs/orientagraph.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["runtime"]
    shell:
        """
        orientagraph \
            -i {input.frq} \
            -o {params.prefix} \
            -m {wildcards.m} \
            -k {params.block_size} \
            {params.root_flag} \
            {params.global_flag} \
            -bootstrap \
            {params.mlno_flag} \
            {params.allmigs_flag} \
            {params.seed_flag} \
            > {log} 2>&1
        """


rule treemix_migration_summary:
    """
    Collect OrientAGraph likelihood files across migration-edge counts.
    """
    input:
        lambda wildcards: expand(
            "results/{project}/treemix/{project}.treemix.m{m}.llik",
            project=wildcards.project,
            m=_treemix_migration_edges(wildcards),
        )
    output:
        summary="results/{project}/treemix/{project}.treemix.migration_summary.tsv"
    run:
        import os

        os.makedirs(os.path.dirname(output.summary), exist_ok=True)
        with open(output.summary, "w") as out:
            out.write("migration_edges\tllik_file\tllik\n")
            for path in input:
                stem = os.path.basename(path)
                m_value = stem.split(".treemix.m", 1)[1].split(".llik", 1)[0]
                with open(path) as handle:
                    llik = " ".join(line.strip() for line in handle if line.strip())
                out.write(f"{m_value}\t{path}\t{llik}\n")


rule treemix_bootstrap_summary:
    """
    Summarize selected-m OrientAGraph / TreeMix bootstrap replicates.
    """
    input:
        llik=lambda wildcards: expand(
            "results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.llik",
            project=wildcards.project,
            m=_treemix_bootstrap_migration_edges(wildcards),
            rep=_treemix_bootstrap_replicates(wildcards),
        ),
        edges=lambda wildcards: expand(
            "results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.edges.gz",
            project=wildcards.project,
            m=_treemix_bootstrap_migration_edges(wildcards),
            rep=_treemix_bootstrap_replicates(wildcards),
        ),
        vertices=lambda wildcards: expand(
            "results/{project}/treemix/bootstrap/{project}.treemix.bootstrap.r{rep}.m{m}.vertices.gz",
            project=wildcards.project,
            m=_treemix_bootstrap_migration_edges(wildcards),
            rep=_treemix_bootstrap_replicates(wildcards),
        )
    output:
        summary="results/{project}/treemix/{project}.treemix.bootstrap_summary.tsv",
        migration_edges="results/{project}/treemix/{project}.treemix.bootstrap_migration_edges.tsv",
        pdf="results/{project}/treemix/plots/{project}.treemix.bootstrap_likelihoods.pdf",
        rds="results/{project}/treemix/plots/{project}.treemix.bootstrap_likelihoods.rds"
    params:
        width=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("likelihood_width", 7),
        height=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("likelihood_height", 5),
        dpi=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("dpi", 300)
    log:
        "logs/{project}/treemix_bootstrap_summary.log"
    benchmark:
        "benchmarks/{project}/treemix_bootstrap_summary.txt"
    conda:
        "../envs/r-ggtree.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/summarize_treemix_bootstrap.R"


rule plot_treemix:
    """
    Plot an OrientAGraph / TreeMix-compatible population graph and residual covariance heatmap.
    """
    input:
        treeout=rules.treemix.output.treeout,
        vertices=rules.treemix.output.vertices,
        edges=rules.treemix.output.edges,
        cov=rules.treemix.output.cov,
        modelcov=rules.treemix.output.modelcov,
        covse=rules.treemix.output.covse,
        populations=rules.treemix_prepare_input.output.populations
    output:
        graph_pdf="results/{project}/treemix/plots/{project}.treemix.m{m}.graph.pdf",
        graph_rds="results/{project}/treemix/plots/{project}.treemix.m{m}.graph.rds",
        residual_pdf="results/{project}/treemix/plots/{project}.treemix.m{m}.residuals.pdf",
        residual_rds="results/{project}/treemix/plots/{project}.treemix.m{m}.residuals.rds"
    params:
        m=lambda wildcards: wildcards.m,
        width=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("width", 10),
        height=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("height", 7),
        dpi=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("dpi", 300),
        migration_arrow_length=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("migration_arrow_length", 0.15)
    log:
        "logs/{project}/plot_treemix.m{m}.log"
    benchmark:
        "benchmarks/{project}/plot_treemix.m{m}.txt"
    conda:
        "../envs/r-ggtree.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_treemix.R"


rule plot_treemix_likelihood:
    """
    Plot final log likelihood against migration-edge count.
    """
    input:
        summary=rules.treemix_migration_summary.output.summary
    output:
        pdf="results/{project}/treemix/plots/{project}.treemix.likelihoods.pdf",
        rds="results/{project}/treemix/plots/{project}.treemix.likelihoods.rds"
    params:
        width=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("likelihood_width", 7),
        height=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("likelihood_height", 5),
        dpi=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("dpi", 300)
    log:
        "logs/{project}/plot_treemix_likelihood.log"
    benchmark:
        "benchmarks/{project}/plot_treemix_likelihood.txt"
    conda:
        "../envs/r-ggtree.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_treemix_likelihood.R"


rule plot_treemix_optm:
    """
    Run OptM on OrientAGraph / TreeMix-compatible outputs and plot migration-edge support.
    """
    input:
        summary=rules.treemix_migration_summary.output.summary,
        llik=lambda wildcards: expand(
            "results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.llik",
            project=wildcards.project,
            m=_treemix_migration_edges(wildcards),
            rep=_treemix_optm_replicates(wildcards),
        ),
        cov=lambda wildcards: expand(
            "results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.cov.gz",
            project=wildcards.project,
            m=_treemix_migration_edges(wildcards),
            rep=_treemix_optm_replicates(wildcards),
        ),
        modelcov=lambda wildcards: expand(
            "results/{project}/treemix/optm_runs/{project}.treemix.r{rep}.m{m}.modelcov.gz",
            project=wildcards.project,
            m=_treemix_migration_edges(wildcards),
            rep=_treemix_optm_replicates(wildcards),
        )
    output:
        tsv="results/{project}/treemix/{project}.treemix.optm_evanno.tsv",
        pdf="results/{project}/treemix/plots/{project}.treemix.optm_evanno.pdf",
        rds="results/{project}/treemix/plots/{project}.treemix.optm_evanno.rds"
    params:
        folder="results/{project}/treemix/optm_runs",
        method=lambda wildcards: _treemix_params(wildcards).get("optm", {}).get("method", "Evanno"),
        width=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("likelihood_width", 7),
        height=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("likelihood_height", 5),
        dpi=lambda wildcards: _treemix_params(wildcards).get("plot", {}).get("dpi", 300)
    log:
        "logs/{project}/plot_treemix_optm.log"
    benchmark:
        "benchmarks/{project}/plot_treemix_optm.txt"
    conda:
        "../envs/r-ggtree.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/plot_treemix_optm.R"
