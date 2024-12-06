# ATAC-seq Analysis Pipeline

## Overview
This pipeline processes ATAC-seq data, including trimming, mapping, peak calling, and quality control. It integrates CHIP-Atlas data and generates final outputs like peaks and normalized BigWig files.

---

## Requirements

### Software:
- **cutadapt** - For adapter trimming.
- **bwa** - For mapping reads to the reference genome.
- **samtools** - For BAM processing.
- **macs2** - For peak calling.
- **aria2c** - For threaded downloads from CHIP-Atlas.
- **deeptools** - For calculating normalized BigWig files and FRiP scores.
- **fastqc** - For FASTQ quality control.
- **multiqc** - For aggregating QC reports.

### Python Libraries:
- `pandas`
- `pysam`
- `deeptools`

---

## Folder Structure

- **config/**: Configuration files (e.g., `config.yaml`).
- **data/**: Raw FASTQ files or SRR files for download.
- **trimmed/**: Trimmed FASTQ files after adapter removal.
- **results/**: Processed outputs including mapping, peaks, and bigwig files.
- **qc/**: Quality control outputs (e.g., FASTQC, flagstat, stats, and FRiP).

---

## Input File
The pipeline requires a **sample metadata file** in TSV format. Below is an example:

| Name       | Unit | Library | Fastq1                                        | Fastq2                                        | GSM |
|------------|------|---------|-----------------------------------------------|-----------------------------------------------|-----|
| LNCaP_ATAC | r1   | Paired  | /groups/lackgrp/raw_data/LNCaP/ATAC-seq/ATAC_16h_R1_1.fq.gz | /groups/lackgrp/raw_data/LNCaP/ATAC-seq/ATAC_16h_R1_2.fq.gz | -   |
| LNCaP_ATAC | r2   | Paired  | /groups/lackgrp/raw_data/LNCaP/ATAC-seq/ATAC_16h_R2_1.fq.gz | /groups/lackgrp/raw_data/LNCaP/ATAC-seq/ATAC_16h_R2_2.fq.gz | -   |

- **Name**: Sample name.
- **Unit**: Unit or replicate identifier.
- **Library**: Single or paired-end sequencing.
- **Fastq1**: Path to raw FASTQ file (R1).
- **Fastq2**: Path to raw FASTQ file (R2 for paired-end data, `-` for single-end).
- **GSM**: GEO sample ID (optional, used for CHIP-Atlas integration).

---

## Configurations
The configuration file (`config/config.yaml`) specifies pipeline parameters. Below is an example:

```yaml
SAMPLES: config/samples.tsv

OUTPUT:
  REF: hg38
  RUN:
    QC: True
    PEAKS: True
    BWS: True
    CHIPATLASBED: False
    CHIPATLASBIGWIG: False
  BW_NORMALIZATIONS:
    - rawcount
    - FPM
  BAMPROCESS_PARAMS: -q 30
  MACS_THRESHOLD: 0.01
  IDR_THRESHOLD: 0.05
  CHIPATLASBED_THRESHOLD: '05'

REFERENCES:
  hg38:
    FA: /groups/lackgrp/genomeAnnotations/hg38/hg38.fa
    BWA_IDX: /groups/lackgrp/genomeAnnotations/hg38/hg38.bwa.idx
    CHROM_SIZES: /groups/lackgrp/genomeAnnotations/hg38/hg38.chrom.sizes
    BLACKLIST: /groups/lackgrp/genomeAnnotations/hg38/hg38-blacklist.v2.bed

```

---

## Running the Pipeline

### 1. Configure the Pipeline
Edit the `config/config.yaml` file to specify:
- Reference genome details.
- Output directory structure.
- Flags to enable or disable specific steps (e.g., QC, peak calling, BigWig normalization).

Ensure that the `config/samples.tsv` file is properly formatted with the sample information.

---

### 2. Slurm Profile
The pipeline uses a Slurm cluster via the `profile/` directory. The `config.yaml` for the Slurm profile should include the following:

```yaml
cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --partition={resources.partition}
    --cpus-per-task={threads}
    --mem={resources.mem_mb}
    --job-name=smk-{rule}-{wildcards}
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
default-resources:
  - partition=normal,big-mem,long,express
  - mem_mb=700000
  - disk_mb=1024000
restart-times: 1
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 12
keep-going: True
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
use-conda: True

```
- **Cluster Resources**: Adjust memory (mem_mb), disk (disk_mb), and partition names according to your Slurm setup.
- **Logging**: Logs for each rule are stored in logs/{rule}/.

---
### 3. Submission Script

Use the following run_pipeline.sh script to submit the pipeline to the Slurm cluster. The script activates the required conda environment and runs Snakemake with the specified profile.

```bash
#!/bin/bash
#SBATCH -c 64
#SBATCH --mem 720GB
#SBATCH -p long,big-mem,normal,express

source ~/.bashrc
conda activate atacseq

snakemake --profile profile/

```
---
### 4. Submit the Pipeline

Run the following command to execute the pipeline:

```bash
sbatch run_pipeline.sh
```
This will:
- Automatically submit jobs to the Slurm cluster.
- Use the configuration specified in the profile/config.yaml file.
- Execute all defined rules in the pipeline.