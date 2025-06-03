rule map_bwa:
    input:
        get_bwa
    output:
        raw=temp("results_{ref}/mapping/{raw}.raw.bam"),
        log="qc/bwa/{ref}:{raw}.bwa.log"

    params:
        idx=lambda wildcards: references[wildcards.ref]["BWA_IDX"]
    threads: 
        32
    conda:
        "../envs/atac.yaml"
    shell:
        """
        bwa mem -t {threads} {params.idx} {input} 2> {output.log} \
        | samtools view -bS - > {output.raw}
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
        "../envs/atac.yaml"
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
        "../envs/atac.yaml"
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
rule bam_merge:
    input:
        get_units          # list of replicate BAMs
    output:
        "results_{ref}/mapping/{name}.final.bam"
    threads: 32
    conda: "../envs/atac.yaml"
    shell:
        r"""
        set -euo pipefail
        bam_out={output}
        if [ $(echo {input} | wc -w) -gt 1 ]; then
            # Multiple replicates → merge
            samtools merge -@ {threads} -o $bam_out {input}
        else
            # Single replicate → just rename
            mv {input} $bam_out
        fi
        samtools index $bam_out
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