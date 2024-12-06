rule macs2_call_peaks:
    input:
        "results_{ref}/mapping/{raw}.filtered.bam"
    output:
        "results_{ref}/peaks/{raw}_peaks.narrowPeak",
        "results_{ref}/peaks/{raw}_peaks.summits.bed"
    params:
        q=config["OUTPUT"]["MACS_THRESHOLD"]
    threads: 16
    shell:
        """
        macs2 callpeak \
            -t {input} \
            -n results_{wildcards.ref}/peaks/{wildcards.raw} \
            -g hs -q {params.q} --nomodel --shift -100 --extsize 200 --outdir results_{wildcards.ref}/peaks
        """

        
rule idr_analysis:
    input:
        get_idr_i  # Dynamically fetch IDR inputs
    output:
        "results_{ref}/idr/{raw}_idr.narrowPeak"
    params:
        idr_threshold=config["OUTPUT"]["IDR_THRESHOLD"]
    shell:
        """
        idr --samples {input['pr1']} {input['pr2']} \
            --peak-list {input['final']} \
            --input-file-type narrowPeak \
            --output-file {output} \
            --rank signal.value \
            --soft-idr-threshold {params.idr_threshold} \
            --plot --use-best-multisummit-IDR
        """