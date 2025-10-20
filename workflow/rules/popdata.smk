rule generate_popmap:
    input:
        vcf=rules.sort_vcf.output.vcf,
        popdata=lambda wildcards: config.get("popdata", "")
    output:
        config["analysis_name"] + "popmap.txt"
    log:
        config["analysis_name"] + "/logs/generate_popmap.log"
    benchmark:
        config["analysis_name"] + "/benchmarks/generate_popmap.txt"
    conda:
        "../envs/bcftools.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"],
        runtime = config["resources"]["default"]["runtime"]
    run:
        import os
        import subprocess
        
        popdata_file = config.get("popdata", "")
        separator = config.get("popseparator", "-")
        
        if popdata_file and os.path.exists(popdata_file):
            # Case 1: Use existing popdata file
            shell("cp {popdata_file} {output.popmap}")
        else:
            # Case 2: Extract population from sample names in VCF
            # Using bcftools to get sample names
            result = subprocess.run(
                ["bcftools", "query", "-l", input.vcf],
                capture_output=True, text=True, check=True
            )
            samples = result.stdout.strip().split("\n")
            
            with open(output.popmap, "w") as f:
                for s in samples:
                    if separator in s:
                        pop = s.split(separator)[0]
                    else:
                        pop = "UNKNOWN"
                    f.write(f"{s}\t{pop}\n")