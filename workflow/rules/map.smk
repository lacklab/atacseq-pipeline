#rule map_bwa:
#    input:
#        get_bwa
#    output:
#        raw=temp("results_{ref}/mapping/{raw}.raw.bam"),
#        log="qc/bwa/{ref}:{raw}.bwa.log"
#
#    params:
#        idx=lambda wildcards: references[wildcards.ref]["BWA_IDX"]
#    threads: 
#        32
#    conda:
#        "../envs/atac.yaml"
#    shell:
#        """
#        bwa mem -t {threads} {params.idx} {input} 2> {output.log} \
#        | samtools view -bS - > {output.raw}
#        """



rule map_bwa:
    input:
        trimmed_fq1="trimmed/{raw}_1.trimmed.fastq.gz",
        trimmed_fq2="trimmed/{raw}_2.trimmed.fastq.gz"
    output:
        raw=temp("results_{ref}/mapping/{raw}.raw.bam"),  # Temporary raw BAM file
        log="qc/bwa/{ref}:{raw}.bwa.log"
    params:
        idx=lambda wildcards: references[wildcards.ref]["BWA_IDX"],
        lib=lambda wildcards: get_lib(wildcards)
    threads:
        16  # Number of threads to use
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        if [[ "{params.lib}" == "Single" ]]; then
            bwa mem -t {threads} {params.idx} {input.trimmed_fq1} 2> {output.log} \
            | samtools view -bS - > {output.raw}
        elif [[ "{params.lib}" == "Paired" ]]; then
            bwa mem -t {threads} {params.idx} {input} 2> {output.log} \
            | samtools view -bS - > {output.raw}
        fi
        """
    # Rule: Process BAM file (coordinate sorting, fixing mates, marking duplicates)
rule bam_process:
    input:
        "results_{ref}/mapping/{raw}.raw.bam"  # Raw BAM file from mapping
    output:
        temp("results_{ref}/mapping/{raw}.coorsorted.bam")  # Coordinate-sorted BAM file
    params:
        config["OUTPUT"]["BAMPROCESS_PARAMS"]  # Parameters for BAM processing
    threads:
        16
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        samtools view -h {params} {input} \
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
    conda:
        "../envs/bedtools.yaml"
    params:
        fa = lambda wildcards: references[wildcards.ref]["FA"],
        bl = lambda wildcards: references[wildcards.ref]["BLACKLIST"],
        p  = lambda wildcards: "-F 3852 -f 2" if get_lib(wildcards) == "Paired" else "-F 3852"
    shell:
        """
        samtools view {input} | egrep -v "chrM" | \
        samtools view -b -@ {threads} -T {params.fa} {params.p} | \
        bedtools intersect -nonamecheck -v -abam stdin -b {params.bl} > {output}

        samtools index {output}
        """

##############################################################################
# Merge technical replicates (or move a single BAM)
##############################################################################
#rule bam_merge:
#    input:
#        get_units          # list of replicate BAMs
#    output:
#        "results_{ref}/mapping/{name}.final.bam"
#    threads: 32
#    conda: "../envs/bwa.yaml"
#    shell:
#        r"""
#        set -euo pipefail
#        bam_out={output}
#        if [ $(echo {input} | wc -w) -gt 1 ]; then
#            # Multiple replicates → merge
#            samtools merge -@ {threads} -o $bam_out {input}
#        else
#            # Single replicate → just rename
#            mv {input} $bam_out
#        fi
#        samtools index $bam_out
#        """
#

rule bam_merge:
    input:
        get_units
    output:
        "results_{ref}/mapping/{name}.final.bam"  # Final merged BAM file
    threads:
        16
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        if [[ $(echo {input} | grep ' ') ]]; then
            samtools merge -@ {threads} -o {output} {input}
        else
            mv {input} {output}
        fi
        samtools index {output}
        """
#rule pseudoreps:
#    input:
#        "results_{ref}/mapping/{raw}.final.bam"
#    output:
#        pr1="results_{ref}/mapping/{raw}.pr1.bam",
#        pr2="results_{ref}/mapping/{raw}.pr2.bam"
#    threads: 16
#    shell:
#        """
#        samtools index {input}
#        samtools view -b --subsample 0.5 {input} > {output.pr1}
#        samtools view -b --subsample 0.5 {input} > {output.pr2}
#        """