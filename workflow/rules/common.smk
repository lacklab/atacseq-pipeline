import pandas as pd
import urllib.request
import os

import yaml
from pathlib import Path



import yaml
from pathlib import Path
references = yaml.safe_load(Path("config/references.yaml").read_text())

# Load sample metadata
samples = pd.read_table(config["SAMPLES"])
samples["Raw"] = samples["Name"] + "." + samples["Unit"].astype(str)

references = yaml.safe_load(Path("config/references.yaml").read_text())


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
#ref = config["OUTPUT"]["REF"]
#macs_threshold = config["OUTPUT"]["MACS_THRESHOLD"]

# Utility functions
def get_fastq(wildcards, fastq_col):
    return samples.loc[samples["Raw"] == wildcards.raw, fastq_col].unique()[0]


def get_lib(wildcards):
    return samples.loc[samples["Raw"] == wildcards.raw, "Library"].unique()[0]

def get_units(wildcards):
    return expand("results_{{ref}}/mapping/{rep}.filtered.bam", rep=samples.loc[samples["Name"] == wildcards.name, "Raw"].unique()) 


def get_group_deeptools(wildcards):
    return expand("results_{{ref}}/mapping/{name}.filtered.bam", name=samples["Raw"].tolist())




# map.smk functions
def get_bwa(wildcards):
    lib = get_lib(wildcards)
    if lib == 'Single':
        return f"trimmed/{wildcards.raw}_1.trimmed.fastq.gz"
    elif lib == 'Paired':
        return (f"trimmed/{wildcards.raw}_1.trimmed.fastq.gz", f"trimmed/{wildcards.raw}_2.trimmed.fastq.gz")
    fq1 = get_fastq(wildcards, "Fastq1")
    
def get_fqs(wildcards):
    fq1 = get_fastq(wildcards, "Fastq1")
    lib = get_lib(wildcards)
    if "SRR" in fq1:
        return (f"sra-data/{fq1}_1.fastq.gz", f"sra-data/{fq1}_2.fastq.gz") if lib == "Paired" else f"sra-data/{fq1}_1.fastq.gz"
    fq2 = get_fastq(wildcards, "Fastq2")
    return (fq1, fq2) if lib == "Paired" else fq1



##########################



# Function to collect outputs for MultiQC
def get_multiqc(wildcards):
    out = []
    
    # List of common quality control file types and tools
    qc_tools_1 = {
        "fastqc": [
            "{raw}_1_fastqc.html", 
            "{raw}_2_fastqc.html"
        ],
        "trimgalore": [
            "{raw}_1.fastq.gz_trimming_report.txt",
            "{raw}_2.fastq.gz_trimming_report.txt",
            "{raw}_1.trimmed_fastqc.html",
            "{raw}_2.trimmed_fastqc.html"
        ],
        "samtools": [
            "flagstat/{ref}:{raw}.coorsorted.flagstat",
            "idxstats/{ref}:{raw}.coorsorted.idxstats",
            "stats/{ref}:{raw}.coorsorted.stats"
        ],
        "bwa": [
            "{ref}:{raw}.bwa.log"
        ],
        "deeptools": [
            "{ref}:all_bam.bamSummary.npz",
            "{ref}:all_bam.plotCorrelation.mat.tab",
            "{ref}:all_bam.plotFingerprint.qcmetrics.txt",
            "{ref}:all_bam.plotFingerprint.raw.txt",
            "{ref}:all_bam.plotPCA.tab"
        ]
    }
    qc_tools_2 = {
         "samtools": [
            "flagstat/{ref}:{name}.final.flagstat",
            "idxstats/{ref}:{name}.final.idxstats",
            "stats/{ref}:{name}.final.stats"
        ],
        "macs": [
            "{ref}:{name}_{q}_peaks.xls"
        ],
        "homer":
        [
            "{ref}:{name}_{q}_summary_mqc.txt"
        ] # TODO: duplications
    }
    # Iterate through each sample and append all files based on the defined templates
    for _, row in samples.iterrows():
        for q in config['OUTPUT']['MACS_THRESHOLD']:
            ref = row['Genome']
            
            raw = row['Raw']
            name = row['Name']

            
            # Generate output paths for each tool and file pattern
            for tool, patterns in qc_tools_1.items():
                for pattern in patterns:
                    out.append(f"qc/{tool}/{pattern.format(raw=raw, ref=ref, q=q)}")

            for tool, patterns in qc_tools_2.items():
                for pattern in patterns:
                    out.append(f"qc/{tool}/{pattern.format(name=name, ref=ref, q=q)}")
   
    # Add FRIP score file (outside the loop as a single file)
        out.append(f"qc/{ref}:frip_mqc.tsv")

    return expand(out)



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
        f"results_{row['Genome']}/peaks/{row['Name']}_{q}_peaks.narrowPeak"
        for i, row in samples.iterrows()
        for q in config['OUTPUT']['MACS_THRESHOLD']
    ]

# BigWig outputsI 
if config["OUTPUT"]["RUN"]["BWS"]:
    outputs += [
        f"results_{row['Genome']}/bigwig/{row['Name']}.bamCoverage.{norm}.bw"
        for i, row in samples.iterrows()
        for norm in config["OUTPUT"]["BW_NORMALIZATIONS"]
    ]