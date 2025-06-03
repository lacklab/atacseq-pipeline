# Rule: Generate alignment statistics using samtools stats
rule stats:
    input:
        "results_{ref}/mapping/{raw}.bam"  # Input BAM file
    output:
        "qc/samtools/stats/{ref}:{raw}.stats"  # Output stats file
    params:
        fa = lambda wildcards: references[wildcards.ref]["FA"]  # Reference genome path
    conda:
        "../envs/atac.yaml"
    shell:
        """
        samtools stats --reference {params.fa} {input} > {output}
        """

# Rule: Generate flagstat statistics using samtools flagstat
rule flagstat:
    input:
        "results_{ref}/mapping/{raw}.bam"  # Input BAM file
    output:
        "qc/samtools/flagstat/{ref}:{raw}.flagstat"  # Output flagstat file
    conda:
        "../envs/atac.yaml"
    shell:
        """
        samtools flagstat {input} > {output}
        """

# Rule: Generate idxstats statistics using samtools idxstats
rule idxstats:
    input:
        "results_{ref}/mapping/{raw}.bam"  # Input BAM file
    output:
        "qc/samtools/idxstats/{ref}:{raw}.idxstats"  # Output idxstats file
    conda:
        "../envs/atac.yaml"
    shell:
        """
        samtools idxstats {input} > {output}
        """

# Rule: Move MACS2 peak calling QC output
rule macs_qc:
    input:
        "results_{ref}/peaks/{name}_{q}_peaks.xls"  # Input MACS2 peaks file
    output:
        "qc/macs/{ref}:{name}_{q}_peaks.xls"  # QC directory for MACS2 peaks file
    conda:
        "../envs/atac.yaml"
    shell:
        """
        ln -s $(readlink -f {input}) $(readlink -f {output})
        """


##############################################################################
# HOMER annotatePeaks summary for MultiQC
##############################################################################
rule annotatepeaks_qc:
    input:
        "results_{ref}/annot/{name}_{q}_annotatepeaks.txt"
    output:
        "qc/homer/{ref}:{name}_{q}_summary_mqc.txt"
    conda: "../envs/atac.yaml"
    shell:
        r"""
        python - <<'PY'
        import pandas as pd
        from collections import Counter
        from pathlib import Path

        infile  = Path({input[0]!r})
        outfile = Path({output[0]!r})
        outfile.parent.mkdir(parents=True, exist_ok=True)

        header_order = [
            "INTERGENIC", "INTRON ", "PROMOTER-TSS ", "EXON ",
            "3' UTR ", "5' UTR ", "TTS ", "NON-CODING "
        ]

        with outfile.open("w") as f:
            f.write(assets["annotatepeaks"])        # MultiQC header

            if infile.stat().st_size == 0:
                counts = {{k: 0 for k in header_order}}
            else:
                df = pd.read_table(infile, comment="#")
                df["shortAnn"] = (
                    df["Annotation"].str.split("(", 1, expand=True)[0].str.upper()
                )
                counts = Counter(df["shortAnn"])

            for key in header_order:
                f.write(f"{{key}}\t{{counts.get(key, 0)}}\n")
        PY
        """

##############################################################################
# FRiP (Fraction of Reads in Peaks)
##############################################################################
rule frip:
    input:
        bams = expand(
            "results_{{ref}}/mapping/{name}.final.bam",
            name=samples["Name"].tolist()
        ),
        peak = expand(
            "results_{{ref}}/peaks/{name}_{q}_peaks.narrowPeak",
            name=samples["Name"].tolist(),
            q=config['OUTPUT']['MACS_THRESHOLD']
        )
    output:
        "qc/{ref}:frip_mqc.tsv"
    conda: "../envs/atac.yaml"
    shell:
        r"""
        python - <<'PY'
        import numpy as np
        import pysam
        from deeptools.countReadsPerBin import CountReadsPerBin
        from pathlib import Path

        bams   = {input.bams!r}
        peaks  = {input.peak!r}
        out_fp = Path({output[0]!r})

        out_fp.parent.mkdir(parents=True, exist_ok=True)
        with out_fp.open("w") as f:
            # MultiQC header
            f.write("# plot_type: 'generalstats'\n")
            f.write("Sample Name\tFRiP\tNumber of Peaks\tMedian Fragment Length\n")

            for bam_file, peak_file in zip(bams, peaks):
                # -------- FRiP ------------
                cr   = CountReadsPerBin([bam_file], bedFile=[peak_file],
                                        numberOfProcessors=10)
                reads_at_peaks = cr.run().sum()

                bam = pysam.AlignmentFile(bam_file)
                total_reads = bam.mapped
                frip = reads_at_peaks / total_reads

                # -------- peaks count -----
                num_peaks = sum(1 for _ in open(peak_file))

                # -------- fragment length -
                frag_lengths = [
                    abs(r.template_length)
                    for r in bam.fetch(until_eof=True) if r.is_proper_pair
                ]
                med_frag = float(np.median(frag_lengths)) if frag_lengths else 0

                sample_name = Path(peak_file).name.split("_peaks")[0]
                f.write(f"{{sample_name}}\t{{frip:.4f}}\t{{num_peaks}}\t{{med_frag:.2f}}\n")
        PY
        """

# Rule: Generate a MultiQC report
rule multiqc:
    input:
        get_multiqc  # Collect all files for MultiQC
    output:
        "qc/multiqc_report.html"  # MultiQC output report
    conda:
        "../envs/atac.yaml"
    shell:
        """
        cd qc/ && multiqc .
        """

# TODO:
# - Integrate FRiP values into 'multiqc_data/multiqc_general_stats.txt'.
# - Add custom content to the MultiQC report (https://multiqc.info/docs/#custom-content).
# - Consider output examples like those from nf-core's ChIP-seq pipeline.