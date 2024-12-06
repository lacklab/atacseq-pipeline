rule homer_annotatepeaks:
	input:
		"results_{ref}/peaks/{raw}_{q}_peaks.narrowPeak"
	output:
		"results_{ref}/annot/{raw}_{q}.annotatepeaks.txt"
	shell:
		"""
		annotatePeaks.pl {input} {ref} > {output}
		"""
