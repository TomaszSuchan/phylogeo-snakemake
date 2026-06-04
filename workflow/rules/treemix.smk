# TreeMix historical population graph inference


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


rule treemix_prepare_input:
    """
    Convert the final analysis VCF to TreeMix allele-count input.
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
        "../envs/treemix.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/vcf_to_treemix.py"


rule treemix:
    """
    Fit a TreeMix graph with a configured number of migration edges.
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
        bootstrap_flag=lambda wildcards: "-bootstrap" if _treemix_params(wildcards).get("bootstrap", False) else "",
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
        "../envs/treemix.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["treemix"]["runtime"]
    shell:
        """
        treemix \
            -i {input.frq} \
            -o {params.prefix} \
            -m {wildcards.m} \
            -k {params.block_size} \
            {params.root_flag} \
            {params.global_flag} \
            {params.bootstrap_flag} \
            {params.seed_flag} \
            > {log} 2>&1
        """


rule treemix_migration_summary:
    """
    Collect TreeMix likelihood files across migration-edge counts.
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
