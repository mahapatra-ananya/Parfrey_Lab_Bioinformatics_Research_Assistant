#!/bin/bash
set -euxo pipefail

####### 00. Setting Up #######
RAW_SEQS_DIR="/home/parfreylab/Desktop/lab_member_files/KMDavis/FTS_metagenomics/raw_data"
PROJECT_DIR="/home/am404/fucus_metagenomics"
RESULTS_DIR="${PROJECT_DIR}/results"
INITIAL_QC_DIR="${RESULTS_DIR}/initial_qc_results"
TRIMMED_DIR="${RESULTS_DIR}/trimmed_sequences"
TRIMMED_QC_DIR="${RESULTS_DIR}/trimmed_qc_results"
ALIGNED_DIR="${RESULTS_DIR}/aligned_sequences"
NONHOST_DIR="${RESULTS_DIR}/nonhost_sequences"
HOST_REFERENCE_DIR="${PROJECT_DIR}/fucus_ncbi_reference/data/GCA_964200245.2"

mkdir -p "$INITIAL_QC_DIR" "$TRIMMED_DIR" "$TRIMMED_QC_DIR" "$ALIGNED_DIR" "$NONHOST_DIR" 

#conda create -n fucus_metagenomic_env -c conda-forge python=3.11 -y
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda config --set channel_priority strict
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate fucus_metagenomic_env
#mamba install -c conda-forge -c bioconda fastqc -y # version 0.12.1
#mamba install -c conda-forge -c bioconda multiqc -y # version 1.35
#mamba install -c conda-forge -c bioconda trimmomatic -y # version 0.41
#mamba install -c bioconda bowtie2 -y # version 2.5.5
#mamba install -c bioconda -c conda-forge samtools # version 1.24

####### 01. Quality Control of Raw Sequences #######
#fastqc -t 8 "$RAW_SEQS_DIR"/*.fastq.gz -o "$INITIAL_QC_DIR"
#cd "$RESULTS_DIR"/initial_qc_results
#multiqc .

####### 02. Trimming Adapters and Low-Quality Bases #######
ADAPTERS="$(find "$CONDA_PREFIX" -name "NexteraPE-PE.fa" -print -quit)"

for r1 in "$RAW_SEQS_DIR"/*_R1_001.fastq.gz; do
    r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    prefix="$(basename "${r1/_R1_001.fastq.gz/}")"

    echo "Trimming sample: $prefix"

    trimmomatic PE -threads 8 -phred33 \
        "$r1" "$r2" \
        "${TRIMMED_DIR}/${prefix}_R1_paired.fastq.gz" \
        "${TRIMMED_DIR}/${prefix}_R1_unpaired.fastq.gz" \
        "${TRIMMED_DIR}/${prefix}_R2_paired.fastq.gz" \
        "${TRIMMED_DIR}/${prefix}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10:2 \
        LEADING:3 \
        TRAILING:3 \
        HEADCROP:15 \
        TAILCROP:5 \
        SLIDINGWINDOW:4:15 \
        MINLEN:100
done

####### 03. Quality Control of Trimmed Sequences #######
fastqc -t 8 "$TRIMMED_DIR"/*.fastq.gz -o "$TRIMMED_QC_DIR"
cd "$TRIMMED_QC_DIR"
multiqc .

####### 04. Alignment to Host Reference Genome #######
bowtie2-build "$HOST_REFERENCE_DIR"/GCA_964200245.2_PHAEFdis_V1.2_genomic.fna "$HOST_REFERENCE_DIR"/GCA_964200245.2_PHAEFdis_V1.2_genomic_index

for r1 in "$TRIMMED_DIR"/*_R1_paired.fastq.gz; do
    r2="${r1/_R1_paired.fastq.gz/_R2_paired.fastq.gz}"
    prefix="$(basename "${r1/_R1_paired.fastq.gz/}")"
    
    echo "Aligning reads to host reference for sample: $prefix"
    
    bowtie2 -x "${HOST_REFERENCE_DIR}/GCA_964200245.2_PHAEFdis_V1.2_genomic_index" \
    -1 "$r1" \
    -2 "$r2" \
    --un-conc-gz "${NONHOST_DIR}/${prefix}_nonhost_%.fastq.gz" \
    -p 8 \
    2> "${ALIGNED_DIR}/${prefix}_bowtie2_report.txt" | \
    samtools sort -@ 8 \
    -o "${ALIGNED_DIR}/${prefix}_host_aligned_sorted.bam"
done

####### 05. Co-Assembly of Metagenomes #######
# R1_FILES=""
# R2_FILES=""

# for r1 in "$NONHOST_DIR"/*_nonhost_1.fastq.gz; do
#     r2="${r1/_nonhost_1.fastq.gz/_nonhost_2.fastq.gz}"

#     R1_FILES="${R1_FILES:+$R1_FILES,}$r1"
#     R2_FILES="${R2_FILES:+$R2_FILES,}$r2"
# done

# megahit \
#     -1 "$R1_FILES" \
#     -2 "$R2_FILES" \
#     -o "${ASSEMBLY_DIR}/coassembly" \
#     -t 16