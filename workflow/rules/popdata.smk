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
    run:
        import pandas as pd

        # Read the provided popdata file
        popdata_df = pd.read_csv(input.popdata, sep="\t", header=0)

        # Read the generated popmap file
        popmap_df = pd.read_csv(input.popmap, sep="\t", header=None, names=["Individual", "Population"])

        # Merge to ensure consistency
        merged_df = pd.merge(popmap_df, popdata_df, left_on="Population", right_on="Population", how="left")


        # Select relevant columns: individual, population, latitude, longitude, and any additional columns
        output_columns = ["Individual", "Population"] + list(merged_df.columns[3:])
        final_popdata_df = merged_df[output_columns]
        
        # Rename the 'Individual_x' column back to 'Individual'
        final_popdata_df = final_popdata_df.rename(columns={"Individual_x": "Individual"})

        # Save the final popdata file
        final_popdata_df.to_csv(output.popdata, sep="\t", header=True, index=False)