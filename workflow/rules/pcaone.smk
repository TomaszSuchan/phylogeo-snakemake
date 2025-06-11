rule NAME:
    input:
        ""
    output:
        ""
    conda:
        "envs/pcaone.yaml"
    shell:
        "pcaone.sh -i {input} -f tsv -o {output}"
