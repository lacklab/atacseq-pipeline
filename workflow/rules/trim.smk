# Rule: Link raw FASTQ files
rule link:
    input:
        get_fqs
    output:
        link_fq1="link/{raw}_1.fastq.gz",
        link_fq2="link/{raw}_2.fastq.gz"
    threads: 2
    run:
        os.makedirs('link', exist_ok=True)
        lib = get_lib(wildcards)
        fq1 = os.path.abspath(input[0])  # Get full absolute path for fq1
        fq2 = os.path.abspath(input[1]) if len(input) > 1 else ""  # Get full absolute path for fq2 (if exists)

        if lib == "Single":
            shell(f"""
                ln -s {fq1} {output.link_fq1}
                touch {output.link_fq2}  # Placeholder for single-end
            """)
        elif lib == "Paired":
            shell(f"""
                ln -s {fq1} {output.link_fq1}
                ln -s {fq2} {output.link_fq2}
            """)

rule trim_adapters:
    input:
        link_fq1 = "link/{raw}_1.fastq.gz",
        link_fq2 = "link/{raw}_2.fastq.gz"
    output:
        trimmed_fq1 = "trimmed/{raw}_1.trimmed.fastq.gz",
        trimmed_fq2 = "trimmed/{raw}_2.trimmed.fastq.gz",
        fastqc1     = "qc/fastqc/{raw}_1_fastqc.html",
        fastqc2     = "qc/fastqc/{raw}_2_fastqc.html",
        t1fastqc    = "qc/trimgalore/{raw}_1.trimmed_fastqc.html",
        t2fastqc    = "qc/trimgalore/{raw}_2.trimmed_fastqc.html",
        t1fastqc_z  = "qc/trimgalore/{raw}_1.trimmed_fastqc.zip",
        t2fastqc_z  = "qc/trimgalore/{raw}_2.trimmed_fastqc.zip",
        t1report    = "qc/trimgalore/{raw}_1.fastq.gz_trimming_report.txt",
        t2report    = "qc/trimgalore/{raw}_2.fastq.gz_trimming_report.txt"
    threads: 8
    conda: "../envs/trim.yaml"
    params:
        lib = lambda wc: get_lib(wc)          # returns "Single" or "Paired"
    shell:
        r"""
        set -euo pipefail

        if [ "{params.lib}" = "Single" ]; then
            ##################################################################
            # SINGLE-END
            ##################################################################
            fastqc --quiet --threads {threads} {input.link_fq1} -o qc/fastqc/

            trim_galore --fastqc --cores {threads} --gzip {input.link_fq1}

            mv {wildcards.raw}_1_trimmed.fq.gz                 {output.trimmed_fq1}
            mv {wildcards.raw}_1_trimmed_fastqc.html           {output.t1fastqc}
            mv {wildcards.raw}_1_trimmed_fastqc.zip            {output.t1fastqc_z}
            mv {wildcards.raw}_1.fastq.gz_trimming_report.txt  {output.t1report}

            # placeholder outputs that don’t exist for SE data
            touch {output.trimmed_fq2} {output.fastqc2} \
                  {output.t2fastqc} {output.t2fastqc_z} {output.t2report}

        else
            ##################################################################
            # PAIRED-END
            ##################################################################
            fastqc --quiet --threads {threads} \
                   {input.link_fq1} {input.link_fq2} -o qc/fastqc/

            trim_galore --fastqc --cores {threads} --paired --gzip \
                        {input.link_fq1} {input.link_fq2}

            mv {wildcards.raw}_1_val_1.fq.gz                   {output.trimmed_fq1}
            mv {wildcards.raw}_2_val_2.fq.gz                   {output.trimmed_fq2}

            mv {wildcards.raw}_1_val_1_fastqc.html             {output.t1fastqc}
            mv {wildcards.raw}_2_val_2_fastqc.html             {output.t2fastqc}

            mv {wildcards.raw}_1_val_1_fastqc.zip              {output.t1fastqc_z}
            mv {wildcards.raw}_2_val_2_fastqc.zip              {output.t2fastqc_z}

            mv {wildcards.raw}_1.fastq.gz_trimming_report.txt  {output.t1report}
            mv {wildcards.raw}_2.fastq.gz_trimming_report.txt  {output.t2report}
        fi
        """