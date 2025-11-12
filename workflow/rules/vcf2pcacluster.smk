# Corrected VCF2PCAcluster rule for multi-project configuration

rule vcf2pcacluster:
    input:
        vcf=rules.thin_vcf.output.vcf  # Use the output from thin_vcf rule
    output:
        eigenvectors = "results/{project}/vcf2pcacluster_miss{miss}_MAF{MAF}/{project}.vcf2pcacluster_miss{miss}_MAF{MAF}.eigenvec",
        eigenvalues = "results/{project}/vcf2pcacluster_miss{miss}_MAF{MAF}/{project}.vcf2pcacluster_miss{miss}_MAF{MAF}.eigenval"
    log:
        "logs/{project}/vcf2pcacluster_miss{miss}_MAF{MAF}.log"
    params:
        # Use wildcards.project instead of config['analysis_name']
        output_prefix=lambda wildcards: (
            f"results/{wildcards.project}/vcf2pcacluster_miss{wildcards.miss}_MAF{wildcards.MAF}/"
            f"{wildcards.project}.vcf2pcacluster_miss{wildcards.miss}_MAF{wildcards.MAF}"
        ),
        MAF=lambda wildcards: wildcards.MAF,
        Miss=lambda wildcards: wildcards.miss,
        # Access project-specific parameters
        cluster_method=lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf2pcacluster"].get("cluster_method", "Kmean"),
        Het=lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf2pcacluster"]["SNP_filtering"].get("Het", 1.0),
        HWE=lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf2pcacluster"]["SNP_filtering"].get("HWE", 0),
        Fchr=lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf2pcacluster"]["SNP_filtering"].get("Fchr", ""),
        KinshipMethod=lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf2pcacluster"].get("KinshipMethod", 1),
        PCnum=lambda wildcards: config["projects"][wildcards.project]["parameters"]["vcf2pcacluster"].get("PCnum", 10)
    wildcard_constraints:
        miss=r"\d+\.?\d*",  # Matches decimal numbers
        MAF=r"\d+\.?\d*"     # Matches decimal numbers
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["vcf2pcacluster"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["vcf2pcacluster"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["vcf2pcacluster"]["runtime"]
    shell:
        """
        test -x {VCF2PCACLUSTER} || {{ echo "ERROR: VCF2PCACluster not found at {VCF2PCACLUSTER}" >&2; exit 1; }}
        VCF2PCACluster -InVCF {input.vcf} \
        -OutPut {params.output_prefix} \
        -Threads {threads} \
        -PCnum {params.PCnum} \
        -Miss {params.Miss} \
        -ClusterMethod {params.cluster_method} \
        -MAF {params.MAF} \
        -Het {params.Het} \
        -HWE {params.HWE} \
        -Fchr {params.Fchr} \
        -KinshipMethod {params.KinshipMethod} \
        &> {log}
        """
