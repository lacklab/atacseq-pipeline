
rule homer_configure:
    output:
        "results_{ref}/.homer_{ref}"
    params:
        genome=lambda wildcards: wildcards.ref
    shell:
        """
        perl /groups/lackgrp/tools/homer/configureHomer.pl -install {params.genome}
        touch {output}
        """

rule homer_annotatepeaks:
    input:
        peaks = "results_{ref}/peaks/{name}_{q}_peaks.narrowPeak",
        prereq = "results_{ref}/.homer_{ref}"
    output:
        "results_{ref}/annot/{name}_{q}_annotatepeaks.txt"
    shell:
        """
        /groups/lackgrp/tools/homer/bin/annotatePeaks.pl {input.peaks} {wildcards.ref} > {output}
        """
