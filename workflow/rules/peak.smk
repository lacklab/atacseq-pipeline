#rule macs2_call_peaks:
#    input:
#        bam="results_{ref}/mapping/{name}.final.bam"
#    output:
#        narrowpeak="results_{ref}/peaks/{name}_{q}_peaks.narrowPeak",
#        summits="results_{ref}/peaks/{name}_{q}_peaks.summits.bed",
#        xls="results_{ref}/peaks/{name}_{q}_peaks.xls"
#    params:
#        p = lambda wildcards: "-f BAMPE" if samples.loc[samples["Name"] == wildcards.name, "Library"].unique()[0] == "Paired" else "-f BAM"
#    threads: 
#        8
#    conda:
#        "../envs/atac.yaml"
#    shell:
#        """
#        macs3 callpeak \
#            -t {input} \
#            -n results_{wildcards.ref}/peaks/{wildcards.name}_{wildcards.q} \
#            -g hs -q {wildcards.q} 
#        """
#
#        #--nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR



def get_macs_p(wildcards):
    lib = samples.loc[samples["Name"] == wildcards.name, "Library"].unique()[0]
    params = f"-t results_{wildcards.ref}/mapping/{wildcards.name}.final.bam "
    params += "-f BAMPE" if lib == "Paired" else "-f BAM"
    return params


# Rule: Call peaks using MACS3
rule macs:
    input:
        "results_{ref}/mapping/{name}.final.bam"  # Function to retrieve MACS input files
    output:
        # MACS3 output files: Excel summary and narrowPeak file
        peaks_xls="results_{ref}/peaks/{name}_{q}_peaks.xls",
        narrow_peak="results_{ref}/peaks/{name}_{q}_peaks.narrowPeak"
    threads:
        4  # Number of threads to use
    params:
        get_macs_p  # Function to retrieve MACS parameters
    conda:
        "../envs/atac.yaml"
    shell:
        """
        macs3 callpeak \
            {params} \
            -g hs \
            -n results_{wildcards.ref}/peaks/{wildcards.name}_{wildcards.q} \
            -q {wildcards.q}  # FDR threshold for peak calling
        """
