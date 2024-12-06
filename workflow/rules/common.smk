import pandas as pd
import urllib.request
import os
import pysam


# Load sample metadata
samples = pd.read_table(config["SAMPLES"])
samples["Raw"] = samples["Name"] + "_" + samples["Unit"].astype(str)

# Asset initialization
assets = {
    "annotatepeaks": "",
    "srx2gsm": {},
    "gsm2srx": {},
}

# Load annotatepeaks asset
with open("assets/annotatepeaks.asset", "r") as f:
    assets["annotatepeaks"] = f.read()

# Download or update CHIP-Atlas experiment list
chip_atlas_file = "assets/experimentList.tab"
if not os.path.exists(chip_atlas_file) or config["UPDATE_CHIPATLAS"]:
    urllib.request.urlretrieve(
        "https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab", chip_atlas_file
    )


# Parse CHIP-Atlas experiment list
with open(chip_atlas_file, "r") as f:
    for line in f:
        cols = line.split("\t")
        if len(cols) < 9 or "SRX" not in cols[0] or cols[8] == "-":
            continue
        gsm, srx = cols[8].split(":")[0], cols[0]
        assets["gsm2srx"][gsm] = srx
        assets["srx2gsm"][srx] = gsm

# Reference genome
ref = config["OUTPUT"]["REF"]
macs_threshold = config["OUTPUT"]["MACS_THRESHOLD"]

# Utility functions
def get_lib(wildcards):
    return samples.loc[samples["Raw"] == wildcards.raw, "Library"].unique()[0]

def get_units(wildcards):
    return samples.loc[samples["Name"] == wildcards.raw, "Raw"].unique()

def get_fq1(wildcards):
    return samples.loc[samples["Raw"] == wildcards.raw, "Fastq1"].unique()[0]

def get_fq2(wildcards):
    return samples.loc[samples["Raw"] == wildcards.raw, "Fastq2"].unique()[0]


def get_fqs(wildcards):
    fq1 = get_fq1(wildcards)
    fq2 = get_fq2(wildcards)
    srr = samples.loc[samples["Raw"] == wildcards.raw, "SRR"].unique()[0]

    # If FASTQ files are not specified, use SRR data
    if fq1 == "-" or fq2 == "-":
        return f"sra-data/{srr}_1.fastq.gz", f"sra-data/{srr}_2.fastq.gz"
    return fq1, fq2

def get_filter_p(wildcards):
    lib = get_lib(wildcards)
    return "-F 3852 -f 2" if lib == "Paired" else "-F 3852"

def get_reps(wildcards):
    return expand(f"results_{ref}/mapping/{{rep}}.filtered.bam", rep=get_units(wildcards))

# QC-related functions
def get_fastqc(wildcards):
    fastqc_map = {}
    for fq in samples["Fastq1"].tolist() + samples["Fastq2"].tolist():
        if fq.endswith(("gz", "zip")):
            fastqc_map[os.path.basename(fq).rsplit(".", 2)[0]] = fq
        elif "SRR" in fq:
            fastqc_map[f"{fq}_1"] = f"sra-data/{fq}_1.fastq.gz"
            fastqc_map[f"{fq}_2"] = f"sra-data/{fq}_2.fastq.gz"
        else:
            fastqc_map[os.path.basename(fq).rsplit(".", 1)[0]] = fq
    return fastqc_map[wildcards.raw]

# Peak-calling functions
def get_macs_p(wildcards):
    lib = get_lib(wildcards)
    control = get_control(wildcards)
    param = f"-t results_{wildcards.ref}/mapping/{wildcards.raw}.{wildcards.p}.bam "
    if control != "-":
        param += f"-c results_{wildcards.ref}/mapping/{control}.final.bam "
    param += "-f BAMPE" if lib == "Paired" else "-f BAM"
    return param

def get_macs_i(wildcards):
    control = get_control(wildcards)
    inputs = [f"results_{wildcards.ref}/mapping/{wildcards.raw}.{wildcards.p}.bam"]
    if control != "-":
        inputs.append(f"results_{wildcards.ref}/mapping/{control}.final.bam")
    return inputs

def get_idr_i(wildcards):
    return {
        "pr1": f"results_{wildcards.ref}/peaks/{wildcards.raw}_pr1_{macs_threshold}_peaks.narrowPeak",
        "pr2": f"results_{wildcards.ref}/peaks/{wildcards.raw}_pr2_{macs_threshold}_peaks.narrowPeak",
        "final": f"results_{wildcards.ref}/peaks/{wildcards.raw}_final_{macs_threshold}_peaks.narrowPeak",
    }

# CHIP-Atlas functions
def get_cabeds(wildcards):
    srx = assets["gsm2srx"][wildcards.gsm]
    return f"results_{wildcards.ref}/cabeds/srx/{srx}.{wildcards.threshold}.bed"

def get_cabws(wildcards):
    srx = assets["gsm2srx"][wildcards.gsm]
    return f"results_{wildcards.ref}/cabws/srx/{srx}.bw"

# Output list
outputs = []

# QC outputs
if config["OUTPUT"]["RUN"]["QC"]:
    outputs.append("qc/multiqc_report.html")

# Peak calling outputs
if config["OUTPUT"]["RUN"]["PEAKS"]:
    idr_threshold = config["OUTPUT"]["IDR_THRESHOLD"]
    outputs += [
        f"results_{ref}/idr/{raw}_{idr_threshold}_idr.narrowPeak"
        for raw in samples["Name"]
        if "input" not in raw.lower()
    ]

# BigWig outputs
if config["OUTPUT"]["RUN"]["BWS"]:
    outputs += [
        f"results_{ref}/bigwig/{raw}.genomecov.{norm}.bw"
        for raw in samples["Name"]
        for norm in config["OUTPUT"]["BW_NORMALIZATIONS"]
    ]