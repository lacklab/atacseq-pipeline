rule macs2_call_peaks:
    input:
        bam="results_{ref}/mapping/{name}.final.bam"
    output:
        narrowpeak="results_{ref}/peaks/{name}_{q}_peaks.narrowPeak",
        summits="results_{ref}/peaks/{name}_{q}_peaks.summits.bed",
        xls="results_{ref}/peaks/{name}_{q}_peaks.xls"
    params:
        p = lambda wildcards: "-f BAMPE" if samples.loc[samples["Name"] == wildcards.name, "Library"].unique()[0] == "Paired" else "-f BAM"
    threads: 
        8
    shell:
        """
        macs2 callpeak \
            -t {input} \
            -n results_{wildcards.ref}/peaks/{wildcards.name}_{wildcards.q} \
            -g hs -q {wildcards.q} --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR
        """