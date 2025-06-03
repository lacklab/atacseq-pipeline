rule sra_prefetch:
    output:
        temp("sra-data/{SRA}/{SRA}.sra")
    conda:
        "../envs/atac.yaml"
    shell:
        """
        prefetch -O sra-data {wildcards.SRA}
        """


rule parallel_fastq_dump:
    input:
        "sra-data/{srr}/{srr}.sra"
    output:
        r1 = "sra-data/{srr}_1.fastq.gz",
        r2 = "sra-data/{srr}_2.fastq.gz"
    threads: 64
    conda:  "../envs/atac.yaml"
    params:
        lib = lambda wc: samples.loc[
            samples["Fastq1"] == wc.srr, "Library"
        ].unique()[0]            # "Single" or "Paired"
    shell:
        r"""
        parallel-fastq-dump -t {threads} --split-files --gzip -s {input} -O sra-data
        # touch a dummy mate-2 file if the run is single-end
        if [ "{params.lib}" = "Single" ]; then
            touch {output.r2}
        fi
        """