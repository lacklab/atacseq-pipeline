rule homer_annotatepeaks:
	input:
		"results_{ref}/peaks/{name}_{q}_peaks.narrowPeak"
	output:
		"results_{ref}/annot/{name}_{q}_annotatepeaks.txt"
	shell:
		"""
		annotatePeaks.pl {input} {wildcards.ref} > {output}
		"""
