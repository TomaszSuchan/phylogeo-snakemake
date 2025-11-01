rule generate_popmap:
    input:
        vcf=rules.sort_vcf.output.vcf
    output:
        popmap="{project}/popmap.txt"
    params:
        popmap=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popmap", ""),
        separator=lambda wildcards: config["projects"][wildcards.project]["parameters"].get("popseparator", "-")
    log:
        "{project}/logs/generate_popmap.log"
    benchmark:
        "{project}/benchmarks/generate_popmap.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    shell:
        r"""
        if [ -n "{params.popmap}" ] && [ -f "{params.popmap}" ]; then
            cp {params.popmap} {output.popmap}
        else
            bcftools query -l {input.vcf} | \
            awk -v sep="{params.separator}" '{{ 
                split($0,a,sep); 
                pop = (length(a)>1) ? a[1] : "UNKNOWN"; 
                print $0 "\t" pop 
            }}' > {output.popmap}
        fi
        """

rule generate_popdata:
    input:
        popmap=rules.generate_popmap.output.popmap
    output:
        indpopdata="{project}/indpopdata.txt"
    params:
        popdata=lambda wildcards: config["projects"][wildcards.project]["parameters"]["popdata"],
    log:
        "{project}/logs/generate_popdata.log"
    benchmark:
        "{project}/benchmarks/generate_popdata.txt"
    conda:
        "../envs/pandas.yaml"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["threads"]
    resources:
        mem_mb=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["mem_mb"],
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["resources"]["default"]["runtime"]
    script:
        "../scripts/generate_popdata.py"