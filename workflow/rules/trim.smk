rule trim_adapters:
    input:
        get_fqs
    output:
        trimmed_fq1="trimmed/{raw}_R1.trimmed.fastq.gz",
        trimmed_fq2="trimmed/{raw}_R2.trimmed.fastq.gz"
    threads: 8
    params:
        adapter_fwd=config["CUT_ADAPTERS"]["ADAPTER_FWD"],
        adapter_rev=config["CUT_ADAPTERS"]["ADAPTER_REV"],
        min_length=config["CUT_ADAPTERS"]["MIN_LENGTH"]
    shell:
        """
        cutadapt \
            -a {params.adapter_fwd} \
            -A {params.adapter_rev} \
            -m {params.min_length} \
            -j {threads} \
            -o {output.trimmed_fq1} \
            -p {output.trimmed_fq2} \
            {input[0]} {input[1]}
        """