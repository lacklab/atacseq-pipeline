import pandas as pd
import urllib.request
import os
import pysam


import yaml
from pathlib import Path
references = yaml.safe_load(Path("config/references.yaml").read_text())

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
if not os.path.exists(chip_atlas_file) or config["CHIPATLAS"]['UPDATE']:
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
#ref = config["OUTPUT"]["REF"]
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
    
    # If FASTQ files are not specified, use SRR data
    if  "SRR" in fq1:
        srr = samples.loc[samples["Raw"] == wildcards.raw, "SRR"].unique()[0]
        return f"sra-data/{srr}_1.fastq.gz", f"sra-data/{srr}_2.fastq.gz"
    return fq1, fq2

def get_filter_p(wildcards):
    lib = get_lib(wildcards)
    return "-F 3852 -f 2" if lib == "Paired" else "-F 3852"

def get_reps(wildcards):
    return expand(f"results_{ref}/mapping/{{rep}}.filtered.bam", rep=get_units(wildcards))

# QC-related functions
# Mapping FASTQ files for FASTQC
def create_fastqc_map(samples):
    """
    Create a mapping of FASTQ files for FASTQC based on the sample information.
    """
    fastqc_map = {}
    for fastq in samples["Fastq1"].tolist() + samples["Fastq2"].tolist():
        if "gz" in fastq or "zip" in fastq:
            key = fastq.rsplit(".", 2)[0].split("/")[-1]
            fastqc_map[key] = fastq
        elif "SRR" in fastq:
            fastqc_map[f"{fastq}_1"] = f"sra-data/{fastq}_1.fastq.gz"
            fastqc_map[f"{fastq}_2"] = f"sra-data/{fastq}_2.fastq.gz"
        elif fastq != "-":
            key = fastq.rsplit(".", 1)[0].split("/")[-1]
            fastqc_map[key] = fastq
    return fastqc_map

# Function to get FASTQC paths
def get_fastqc(wildcards):
    return fastqc_map[wildcards.raw]

# Function to get annotate peaks files
def get_annotatepeaks(wildcards):
    """
    Generate paths for annotatePeaks output files.
    """
    return expand(
        [f"qc/{ref}:{row['Name']}.annotatePeaks.txt" for _, row in samples.iterrows()]
    )

# Function to get BAM files for FRiP
def get_frip_b(wildcards):
    """
    Generate paths for final BAM files used in FRiP calculations.
    """
    return expand(
        [f"results_{ref}/mapping/{row['Name']}.final.bam" for _, row in samples.iterrows()]
    )

# Function to get peak files for FRiP
def get_frip_p(wildcards):
    """
    Generate paths for narrowPeak files used in FRiP calculations.
    """
    return expand(
        [f"results_{ref}/peaks/{row['Name']}_{q}_peaks.narrowPeak" for _, row in samples.iterrows()]
    )

# Function to collect outputs for MultiQC
def get_multiqc(wildcards):
    """
    Generate a list of all files required for MultiQC analysis.
    """
    output_files = []
    for _, row in samples.iterrows():
        lib = row["Library"]
        fq1 = row["Fastq1"].split("/")[-1]
        fq2 = row["Fastq2"].split("/")[-1]

        # Adjust file extensions for SRR and compressed files
        if "SRR" in fq1:
            ext1, ext2 = "_1", "_2"
            fq2 = fq1
        else:
            ext1, ext2 = "", ""
        if "gz" in fq1:
            fq1 = fq1.rsplit(".", 2)[0]
            fq2 = fq2.rsplit(".", 2)[0]

        # Add FASTQC outputs
        if lib == "Single":
            output_files.append(f"qc/fastqc/{fq1 + ext1}_fastqc.zip")
        elif lib == "Paired":
            output_files.extend([
                f"qc/fastqc/{fq1 + ext1}_fastqc.zip",
                f"qc/fastqc/{fq2 + ext2}_fastqc.zip"
            ])

        # Add flagstats and stats outputs
        output_files.extend([
            f"qc/flagstats/{ref}:{row['Raw']}.raw",
            f"qc/flagstats/{ref}:{row['Name']}.final",
            f"qc/stats/{ref}:{row['Name']}.final"
        ])

    # Add additional QC outputs
    for _, row in samples.iterrows():
        output_files.extend([
            f"qc/annot/{ref}:{row['Name']}_{q}.summary_mqc.txt",
            f"qc/macs/{ref}:{row['Name']}_{q}_peaks.xls"
        ])

    # Add FRiP QC summary
    output_files.append("qc/frip_mqc.tsv")
    return expand(output_files)
	
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
    outputs += [
        f"results_{ref}/idr/{raw}_idr.narrowPeak"
        for raw in samples["Name"]
    ]

# BigWig outputsI 
if config["OUTPUT"]["RUN"]["BWS"]:
    outputs += [
        f"results_{ref}/bigwig/{raw}_{norm}.bw"
        for raw in samples["Name"]
        for norm in config["OUTPUT"]["BW_NORMALIZATIONS"]
    ]