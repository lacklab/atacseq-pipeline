rule generate_bigwig:
    input:
        "results_{ref}/mapping/{raw}.final.bam"
    output:
        bigwig="results_{ref}/bigwig/{raw}.bamCoverage.{norm}.bw"
    threads: 
        8
    conda:
        "../envs/atac.yaml"
    shell:
        """
        bamCoverage -b {input} -o {output} \
        --binSize 50 \
        --normalizeUsing {wildcards.norm} \
        --extendReads --centerReads
        """



    