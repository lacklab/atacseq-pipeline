rule map_bwa:
    input:
        fq1="trimmed/{raw}_R1.trimmed.fastq.gz",
        fq2="trimmed/{raw}_R2.trimmed.fastq.gz"
    output:
        temp("results_{ref}/mapping/{raw}.raw.bam")
    params:
        idx=lambda wildcards: config["REFERENCES"][wildcards.ref]["BWA_IDX"]
    threads: 32
    shell:
        """
        bwa mem -t {threads} {params.idx} {input.fq1} {input.fq2} \
        | samtools view -bS - > {output}
        """


rule bam_process:
    input:
        "results_{ref}/mapping/{raw}.raw.bam"
    output:
        temp("results_{ref}/mapping/{raw}.coorsorted.bam")
    threads: 32
    params:
        bam_params=config["OUTPUT"]["BAMPROCESS_PARAMS"]
    shell:
        """
        samtools view -h {params.bam_params} {input} \
        | samtools fixmate -m -@ {threads} - - \
        | samtools sort -@ {threads} -m 10G - \
        | samtools markdup -@ {threads} - {output}
        """

rule bam_filter:
    input:
        "results_{ref}/mapping/{raw}.coorsorted.bam"
    output:
        temp("results_{ref}/mapping/{raw}.filtered.bam")
    threads: 32
    params:
        fa=lambda wildcards: config["REFERENCES"][ref]["FA"],
        blacklist=lambda wildcards: config["REFERENCES"][ref]["BLACKLIST"],
        filter_p=get_filter_p
    shell:
        """
        samtools view {input} | egrep -v "chrM" | \
        samtools view -b -@ {threads} -T {params.fa} {params.filter_p} | \
        bedtools intersect -nonamecheck -v -abam stdin -b {params.blacklist} > {output}
        """

rule bam_merge:
    input:
        get_reps  # Collect replicates
    output:
        "results_{ref}/mapping/{raw}.final.bam"
    threads: 32
    run:
        if len(input) > 1:
            shell("""
                samtools merge -@ {threads} -o {output} {input}
                samtools index {output}
            """)
        else:
            shell("""
                mv {input[0]} {output}
                samtools index {output}
            """)

rule pseudoreps:
    input:
        "results_{ref}/mapping/{raw}.final.bam"
    output:
        pr1="results_{ref}/mapping/{raw}.pr1.bam",
        pr2="results_{ref}/mapping/{raw}.pr2.bam"
    threads: 16
    shell:
        """
        samtools index {input}
        samtools view -b --subsample 0.5 {input} > {output.pr1}
        samtools view -b --subsample 0.5 {input} > {output.pr2}
        """