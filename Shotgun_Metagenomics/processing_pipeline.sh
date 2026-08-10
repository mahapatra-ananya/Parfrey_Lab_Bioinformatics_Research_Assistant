####### 00. Setting Up #######
RAW_SEQS_DIR="/home/parfreylab/Desktop/lab_member_files/KMDavis/FTS_metagenomics/raw_data"
PROJECT_DIR="/home/am404/fucus_metagenomics"
RESULTS_DIR="${PROJECT_DIR}/results"

mkdir -p "$RESULTS_DIR"

conda create -n fucus_metagenomic_env -c conda-forge python=3.11 -y
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda activate fucus_metagenomic_env
mamba install -c conda-forge -c bioconda fastqc -y
mamba install -c conda-forge -c bioconda multiqc -y

####### 01. Quality Control of Raw Sequences #######
mkdir "$RESULTS_DIR"/initial_qc_results
fastqc -t 4 "$RAW_SEQS_DIR"/*.fastq.gz -o "$RESULTS_DIR"/qc_results
cd "$RESULTS_DIR"/initial_qc_results
multiqc .


####### 02. Align to Host Reference Genome #######


