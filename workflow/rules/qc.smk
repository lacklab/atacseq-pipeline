# Rule: Generate alignment statistics using samtools stats
rule stats:
    input:
        "results_{ref}/mapping/{raw}.bam"  # Input BAM file
    output:
        "qc/samtools/stats/{ref}:{raw}.stats"  # Output stats file
    params:
        fa = lambda wildcards: references[wildcards.ref]["FA"]  # Reference genome path
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
    shell:
        """
        ln -s $(readlink -f {input}) $(readlink -f {output})
        """

# Rule: Annotate peaks and generate a summary for QC
from collections import Counter
rule annotatepeaks_qc:
    input:
        "results_{ref}/annot/{name}_{q}_annotatepeaks.txt"  # Input annotated peaks file
    output:
        "qc/homer/{ref}:{name}_{q}_summary_mqc.txt"  # Summary QC file
    run:
        header = ["INTERGENIC", "INTRON ", "PROMOTER-TSS ", "EXON ", "3' UTR ", "5' UTR ", "TTS ", "NON-CODING "]
        with open(output[0], "w") as f:
            f.write(assets["annotatepeaks"])  # Pre-written asset for annotatepeaks
            tmp = pd.read_table(input[0])  # Read annotation file
            if tmp.shape[0] == 0:
                nAnnot = dict(zip(header, [0] * len(header)))  # No annotations found
            else:
                tmp["shortAnn"] = tmp["Annotation"].str.split("(", expand=True)[0].str.upper()
                nAnnot = Counter(tmp["shortAnn"])  # Count annotations
            for k in header:
                f.write(f"{k}\t{nAnnot.get(k, 0)}\n")  # Write counts to summary file

# Rule: Compute FRiP (Fraction of Reads in Peaks) values


rule frip:
    input:
        bams=expand("results_{{ref}}/mapping/{name}.final.bam", name=samples["Name"].tolist()),   # Fetch BAM files for FRIP calculation
        peak=expand("results_{{ref}}/peaks/{name}_{q}_peaks.narrowPeak", name=samples["Name"].tolist(), q=config['OUTPUT']['MACS_THRESHOLD'])    # Fetch peak files for FRIP calculation
    output:
        "qc/{ref}:frip_mqc.tsv"
    run:
        import deeptools.countReadsPerBin as crpb
        import pysam
        import numpy as np

        # Prepare the MultiQC-compatible header
        with open(output[0], "w") as f:
            f.write("# plot_type: 'generalstats'\n")
            f.write("Sample Name\tFRiP\tNumber of Peaks\tMedian Fragment Length\n")

            # Loop through BAM and Peak pairs
            for b, p in zip(input.bams, input.peak):
                
                # Calculate FRIP using deepTools
                cr = crpb.CountReadsPerBin([b], bedFile=[p], numberOfProcessors=10)
                reads_at_peaks = cr.run()
                total_reads_at_peaks = reads_at_peaks.sum(axis=0)

                # Calculate total mapped reads using pysam
                bam = pysam.AlignmentFile(b)
                total_mapped_reads = bam.mapped

                # Calculate number of peaks
                with open(p, 'r') as peak_file:
                    num_peaks = sum(1 for _ in peak_file)

                # Calculate median fragment length using pysam
                fragment_lengths = [
                    abs(read.template_length) for read in bam.fetch() if read.is_proper_pair
                ]
                median_fragment_length = np.median(fragment_lengths)

                # Calculate FRIP score
                sample_name = p.split("/")[-1].split("_peaks")[0]
                frip_score = float(total_reads_at_peaks[0]) / total_mapped_reads

                # Write results into a MultiQC-compatible TSV file
                f.write(f"{sample_name}\t{frip_score:.4f}\t{num_peaks}\t{median_fragment_length:.2f}\n")

# Rule: Generate a MultiQC report
rule multiqc:
    input:
        get_multiqc  # Collect all files for MultiQC
    output:
        "qc/multiqc_report.html"  # MultiQC output report
    shell:
        """
        cd qc/ && multiqc .
        """

# TODO:
# - Integrate FRiP values into 'multiqc_data/multiqc_general_stats.txt'.
# - Add custom content to the MultiQC report (https://multiqc.info/docs/#custom-content).
# - Consider output examples like those from nf-core's ChIP-seq pipeline.