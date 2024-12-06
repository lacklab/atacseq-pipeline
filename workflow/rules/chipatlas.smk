rule ChipAtlasBed:
    output:
        "results_{ref}/cabeds/{raw}.{gsm}.{threshold}.bed"
    threads:
        16
    shell:
        """
        url=https://chip-atlas.dbcls.jp/data/{wildcards.ref}/eachData/bed{wildcards.threshold}/{wildcards.srx}.{wildcards.threshold}.bed

        if wget --spider $url 2>/dev/null; then
            /home/ualtintas/apps/aria2-1.35.0/src/aria2c -x {threads} -s {threads} \
                $url \
                -o {output}
        else
            touch {output}  # Create an empty file if the URL doesn't exist
        fi
        """


rule ChipAtlasBigwig:
    output:
        "results_{ref}/cabws/{raw}.{gsm}.bw"
    threads:
        16
    shell:
        """
        url=https://chip-atlas.dbcls.jp/data/{wildcards.ref}/eachData/bw/{wildcards.srx}.bw

        if wget --spider $url 2>/dev/null; then
            /home/ualtintas/apps/aria2-1.35.0/src/aria2c -x {threads} -s {threads} \
                $url \
                -o {output}
        else
            touch {output}  # Create an empty file if the URL doesn't exist
        fi
        """