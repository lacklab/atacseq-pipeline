rule macs2_call_peaks:
    input:
        "results_{ref}/mapping/{raw}.{p}.bam"
    output:
        narrowpeak="results_{ref}/peaks/{raw}_{p}_peaks.narrowPeak",
        summits="results_{ref}/peaks/{raw}_{p}_peaks.summits.bed",
        xls="results_{ref}/peaks/{raw}_{p}_peaks.xls"
    params:
        q=config["OUTPUT"]["MACS_THRESHOLD"],
        outdir="results_{ref}/peaks"
    threads: 16
    shell:
        """
        macs2 callpeak \
            -t {input} \
            -n {params.outdir}/{wildcards.raw}_{wildcards.p} \
            -g hs -q {params.q} --nomodel --shift -100 --extsize 200 --outdir {params.outdir}
        """

rule idr_analysis:
    input:
        pr1=lambda wildcards: f"results_{wildcards.ref}/peaks/{wildcards.raw}_pr1_peaks.narrowPeak",
        pr2=lambda wildcards: f"results_{wildcards.ref}/peaks/{wildcards.raw}_pr2_peaks.narrowPeak",
        final=lambda wildcards: f"results_{wildcards.ref}/peaks/{wildcards.raw}_final_peaks.narrowPeak"
    output:
        "results_{ref}/idr/{raw}_idr.narrowPeak"
    params:
        idr_threshold=config["OUTPUT"]["IDR_THRESHOLD"]
    threads: 8
    shell:
        """
        idr --samples {input.pr1} {input.pr2} \
            --peak-list {input.final} \
            --input-file-type narrowPeak \
            --output-file {output} \
            --rank signal.value \
            --soft-idr-threshold {params.idr_threshold} \
            --plot --use-best-multisummit-IDR
        """