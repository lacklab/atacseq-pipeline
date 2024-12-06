rule fastqc:
    input:
        get_fastqc
    output:
        "qc/fastqc/{raw}_fastqc.zip"
    shell:
        """
        fastqc {input} -o qc/fastqc
        """


rule stats:
    input:
        "results_{ref}/mapping/{raw}.final.bam"
    output:
        "qc/stats/{ref}:{raw}.txt"
    shell:
        """
        samtools stats {input} > {output}
        """

rule flagstat:
    input:
        "results_{ref}/mapping/{raw}.final.bam"
    output:
        "qc/flagstats/{ref}:{raw}.txt"
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule macs_qc:
	input:
		"results_{ref}/peaks/{raw}_{q}_peaks.xls"
	output:
		"qc/macs/{ref}:{raw}_{q}_peaks.xls"
	shell:
		"""
		mv {input} {output}
		"""


from collections import Counter
rule annotatepeaks_qc:
	input:
		"results_{ref}/annot/{raw}_{q}.annotatepeaks.txt"
	output:
		"qc/annot/{ref}:{raw}_{q}.summary_mqc.txt"
	run:
		header = ["INTERGENIC", "INTRON ", "PROMOTER-TSS ", "EXON ", "3' UTR ", "5' UTR ", "TTS ", "NON-CODING "]
		with open(output[0], "w") as f:
			f.write(assets["annotatepeaks"])
			tmp = pd.read_table(input[0])
			if tmp.shape[0] == 0:
				nAnnot = dict(zip(header, [8]*0))
			else:
				tmp["shortAnn"] = tmp["Annotation"].str.split("(", expand=True)[0].str.upper()
				nAnnot = Counter(tmp["shortAnn"])
			for k in header:
				f.write(f"{k}\t{nAnnot[k]}\n")


rule frip:
    input:
        bams=get_frip_b,
        peak=get_frip_p
    output:
        "qc/frip_mqc.tsv"
    run:
        from collections import Counter
        import deeptools.countReadsPerBin as crpb
        with open(output[0], "w") as f:
            f.write("# plot_type: 'generalstats'\nSample Name\tFRiP\n")
            for bam, peak in zip(input.bams, input.peak):
                cr = crpb.CountReadsPerBin([bam], bedFile=[peak], numberOfProcessors=10)
                reads_at_peaks = cr.run()
                total_reads = reads_at_peaks.sum(axis=0)
                mapped_reads = pysam.AlignmentFile(bam).mapped
                name = peak.split("/")[-1].split("_peaks")[0]
                frip = float(total_reads[0]) / mapped_reads
                f.write(f"{name}\t{frip}\n")

# TODO: This can be written in a script.	

rule multiqc:
	input:
		get_multiqc
	output:
		"qc/multiqc_report.html"
	shell:	
		"""
		cd qc/ && multiqc .
		"""

# TODO: get FRiP values and add to 'multiqc_data/multiqc_general_stats.txt'
# TODO: work on custom-content (https://multiqc.info/docs/#custom-content)
# TODO: QC output can include outputs; such as https://nf-co.re/chipseq/1.2.2/output