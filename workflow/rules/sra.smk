rule parallel_fastq_dump:
    output:
        r1 = "sra-data/{srr}_1.fastq.gz",
        r2 = "sra-data/{srr}_2.fastq.gz"
    conda:
        "../envs/sra.yaml"
    threads:
        16
    resources:
        slurm_partition='express,normal,long,big-mem'
    shell:
        r"""
        LIBTYPE=$(awk -v srr={wildcards.srr} '$8==srr{{print $4}}' config/samples.tsv)
        
        echo $LIBTYPE

        if [ "$LIBTYPE" = "Single" ]; then
            parallel-fastq-dump -t {threads} --split-files --gzip -s {wildcards.srr} -O sra-data
            touch {output.r2}  # create dummy r2 for Single-end
        elif [ "$LIBTYPE" = "Paired" ]; then
            parallel-fastq-dump -t {threads} --split-files --gzip -s {wildcards.srr} -O sra-data
        else
            echo "Error: Unknown library type for {wildcards.srr}" >&2
            exit 1
        fi
        """


