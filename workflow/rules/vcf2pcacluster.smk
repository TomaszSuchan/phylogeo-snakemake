# Rule for VCF2PCAcluster
rule vcf2pcacluster:
    input:
        vcf=rules.thin_vcf.output.vcf
    output:
        eigenvectors=config["analysis_name"] + "/vcf2pcacluster_miss{miss}_MAF{MAF}/vcf2pcacluster_miss{miss}_MAF{MAF}.eigenvec",
        eigenvalues=config["analysis_name"] + "/vcf2pcacluster_miss{miss}_MAF{MAF}/vcf2pcacluster_miss{miss}_MAF{MAF}.eigenval"
    log:
        config["analysis_name"] + "/logs/vcf2pcacluster_miss{miss}_MAF{MAF}.log"
    params:
        bin="workflow/bin/VCF2PCACluster",
        output_prefix=lambda wildcards: (
            f"{config['analysis_name']}/vcf2pcacluster_miss{wildcards.miss}_MAF{wildcards.MAF}/"
            f"vcf2pcacluster_miss{wildcards.miss}_MAF{wildcards.MAF}"
        ),
        MAF=lambda wildcards: wildcards.MAF,
        Miss=lambda wildcards: wildcards.miss,
        cluster_method=config["vcf2pcacluster"].get("cluster_method", "Kmean"),
        Het=config["vcf2pcacluster"]["SNP_filtering"].get("Het", 1.0),
        HWE=config["vcf2pcacluster"]["SNP_filtering"].get("HWE", 0),
        Fchr=config["vcf2pcacluster"]["SNP_filtering"].get("Fchr", ""),
        KinshipMethod=config["vcf2pcacluster"].get("KinshipMethod", 1),
        PCnum=config["vcf2pcacluster"].get("PCnum", 10)
    wildcard_constraints:
        miss=r"\d+\.?\d*",  # Matches decimal numbers
        MAF=r"\d+\.?\d*"     # Matches decimal numbers
    threads: config["resources"]["vcf2pcacluster"]["threads"]
    resources:
        mem_mb=config["resources"]["vcf2pcacluster"]["mem_mb"],
        time=config["resources"]["vcf2pcacluster"]["runtime"]
    shell:
        """
        {params.bin} -InVCF {input.vcf} \
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
