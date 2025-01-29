
rule plot_fingerprint:
	input:
		get_group_deeptools
	output:
		plot="qc/deeptools/{ref}:all_bam.plotFingerprint.pdf",
		counts="qc/deeptools/{ref}:all_bam.plotFingerprint.raw.txt",
		qcmetrics="qc/deeptools/{ref}:all_bam.plotFingerprint.qcmetrics.txt"
	threads:
		64
	shell:
		"""	
		plotFingerprint \
			--skipZeros \
			--bamfiles {input} \
			--plotFile {output.plot} \
			--outRawCounts {output.counts} \
			--outQualityMetrics {output.qcmetrics} \
			--numberOfProcessors {threads}
		"""



rule plot_multibam:
	input:
		get_group_deeptools
	output:
		corr="qc/deeptools/{ref}:all_bam.bamSummary.npz",
		plot_pca="qc/deeptools/{ref}:all_bam.plotPCA.pdf",
		tab_pca="qc/deeptools/{ref}:all_bam.plotPCA.tab",
		plot_corr="qc/deeptools/{ref}:all_bam.plotCorrelation.pdf",
		tab_corr="qc/deeptools/{ref}:all_bam.plotCorrelation.mat.tab"
	threads:
		64
	shell:
		"""	
		multiBamSummary bins \
			--smartLabels --binSize 500 \
			--bamfiles {input} \
			--numberOfProcessors {threads} \
			--outFileName {output.corr}

		plotPCA \
			--corData {output.corr} \
			--plotFile {output.plot_pca} \
			--outFileNameData {output.tab_pca}

		plotCorrelation \
			--skipZeros --removeOutliers \
			--corData {output.corr} \
			--corMethod pearson \
			--whatToPlot heatmap \
			--plotFile {output.plot_corr} \
			--outFileCorMatrix  {output.tab_corr}
		"""	