rule generate_popmap:
    input:
        vcf=rules.sort_vcf.output.vcf
    output:
        popmap=config["analysis_name"] + "/popmap.txt"
    params:
        popmap=config.get("popmap", ""),
        separator=config.get("popseparator", "-")
    log:
        config["analysis_name"] + "/logs/generate_popmap.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/generate_popmap.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem_mb"],
        runtime=config["resources"]["default"]["runtime"]
    group: "plots"
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
        popdata=config["popdata"],
        popmap=rules.generate_popmap.output.popmap
    output:
        popdata=config["analysis_name"] + "/popdata.txt"
    log:
        config["analysis_name"] + "/logs/generate_popdata.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/generate_popdata.txt"
    conda:
        "../envs/pandas.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem_mb"],
        runtime=config["resources"]["default"]["runtime"]
    group: "plots"
    script:
        "../scripts/generate_popdata.py"