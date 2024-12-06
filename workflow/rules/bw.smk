rule generate_bigwig:
    input:
        "results_{ref}/mapping/{raw}.final.bam"
    output:
        bigwig="results_{ref}/bigwig/{raw}_{norm}.bw"
    params:
        norm_method="{norm}",
        bin_size=50
    threads: 8
    shell:
        """
        bamCoverage -b {input} -o {output} \
        --binSize {params.bin_size} \
        --normalizeUsing {params.norm_method} \
        --extendReads --centerReads
        """