rule generate_popmap:
    input:
        vcf=rules.sort_vcf.output.vcf
    output:
        popmap=config["analysis_name"] + "/popmap.txt"
    params:
        popdata=config.get("popdata", ""),
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
    shell:
        r"""
        if [ -n "{params.popdata}" ] && [ -f "{params.popdata}" ]; then
            cp {params.popdata} {output.popmap}
        else
            bcftools query -l {input.vcf} | \
            awk -v sep="{params.separator}" '{{ 
                split($0,a,sep); 
                pop = (length(a)>1) ? a[1] : "UNKNOWN"; 
                print $0 "\t" pop 
            }}' > {output.popmap}
        fi
        """