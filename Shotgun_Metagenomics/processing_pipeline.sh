#!/bin/bash
set -euxo pipefail

####### 00. Setting Up #######
RAW_SEQS_DIR="/home/parfreylab/Desktop/lab_member_files/KMDavis/FTS_metagenomics/raw_data"
PROJECT_DIR="/home/am404/fucus_metagenomics"
RESULTS_DIR="${PROJECT_DIR}/results"
INITIAL_QC_DIR="${RESULTS_DIR}/initial_qc_results"
TRIMMED_DIR="${RESULTS_DIR}/trimmed_sequences"
TRIMMED_QC_DIR="${RESULTS_DIR}/trimmed_qc_results"
HOST_ALIGNED_DIR="${RESULTS_DIR}/host_aligned_sequences"
NONHOST_DIR="${RESULTS_DIR}/nonhost_sequences"
HOST_REFERENCE_DIR="${PROJECT_DIR}/fucus_ncbi_reference/data/GCA_964200245.2"
ASSEMBLY_DIR="${RESULTS_DIR}/assembled_sequences"
ALIGNED_DIR="${RESULTS_DIR}/assembly_aligned_sequences"
BIN_DIR="${RESULTS_DIR}/binned_sequences"
TAXONOMY_DIR="${RESULTS_DIR}/taxonomic_classification"

mkdir -p "$INITIAL_QC_DIR" "$TRIMMED_DIR" "$TRIMMED_QC_DIR" "$HOST_ALIGNED_DIR" "$NONHOST_DIR" "$ASSEMBLY_DIR" "$ALIGNED_DIR" "$BIN_DIR" "$TAXONOMY_DIR"

#wget https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz
#tar xzvf master.tar.gz
#wget https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
#tar -xvzf boost_1_86_0.tar.gz
#cd boost_1_86_0/
#./bootstrap.sh --prefix=$HOME/boost-1.86-install
#./b2 install
#export BOOST_ROOT=/home/am404/b2 install
#cd /home/am404/fucus_metagenomics/berkeleylab-metabat-*
#mkdir build && cd build
#cmake \
#  -DBOOST_ROOT=$HOME/boost-1.86-install \
#  -DBoost_NO_SYSTEM_PATHS=ON \
#  -DCMAKE_INSTALL_PREFIX=$HOME/Metabat \
#  ..
#make && make test && make install
#conda create -n fucus_metagenomic_env -c conda-forge python=3.11 -y
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda config --set channel_priority strict
#CONDA_BASE=$(conda info --base)
#source "${CONDA_BASE}/etc/profile.d/conda.sh"
#conda activate fucus_metagenomic_env
#mamba install -c conda-forge -c bioconda fastqc -y # version 0.12.1
#mamba install -c conda-forge -c bioconda multiqc -y # version 1.35
#mamba install -c conda-forge -c bioconda trimmomatic -y # version 0.41
#mamba install -c bioconda bowtie2 -y # version 2.5.5
#mamba install -c bioconda -c conda-forge samtools # version 1.24
#mamba install -c bioconda -c conda-forge megahit # version 1.2.9
#mamba install -c bioconda quast # version 5.3.0
#conda deactivate
#mamba create -n checkm2 -c bioconda -c conda-forge checkm2
#conda activate checkm2
#mkdir -p "${PROJECT_DIR}/checkm2_db"
#checkm2 database --download --path "${PROJECT_DIR}/checkm2_db"
#checkm2 testrun
#conda deactivate
#cd $PROJECT_DIR
#git clone https://github.com/bxlab/metaWRAP.git
#configure metaWRAP/bin/config-metawrap
#echo 'export PATH="$PROJECT_DIR/metaWRAP/bin:$PATH"' >> ~/.bashrc
#mamba create -y -n fucus_metagenomics_metawrap python=2.7
#conda activate fucus_metagenomics_metawrap
#mkdir -p metaWRAP/databases metaWRAP/databases/CHECKM
#cd metaWRAP/databases/CHECKM
#wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
#tar -xvf *.tar.gz
#rm *.gz
#cd ../
# cd $PROJECT_DIR
# mkdir kraken2
# cd $PROJECT_DIR/kraken2
# wget https://github.com/DerrickWood/kraken2/archive/refs/tags/2.17.2.zip
# unzip 2.17.2.zip
# mkdir kraken2_installation
# cd kraken2-2.17.2/
# ./install_kraken2.sh $PROJECT_DIR/kraken2/kraken2_installation
# cp $PROJECT_DIR/kraken2/kraken2_installation/kraken2{,-build,-inspect} $HOME/bin
# #kraken2-build --db $PROJECT_DIR/kraken2/SILVA --special silva --threads 6 # NOT HELPFUL TO USE MARKER GENE
# kraken2-build --standard --threads 6 --db $PROJECT_DIR/kraken2/standard_db --use-ftp
# mamba create -n gtdbtk -c conda-forge -c bioconda gtdbtk=2.7.2
# wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/environment.yaml
# mamba env create -f environment.yaml -n DRAM

CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate fucus_metagenomic_env


####### 01. Quality Control of Raw Sequences #######
fastqc -t 6 "$RAW_SEQS_DIR"/*.fastq.gz -o "$INITIAL_QC_DIR"
cd "$RESULTS_DIR"/initial_qc_results
multiqc .

####### 02. Trimming Adapters and Low-Quality Bases #######
ADAPTERS="$(find "$CONDA_PREFIX" -name "NexteraPE-PE.fa" -print -quit)"

for r1 in "$RAW_SEQS_DIR"/*_R1_001.fastq.gz; do
    r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    prefix="$(basename "${r1/_R1_001.fastq.gz/}")"

    echo "Trimming sample: $prefix"

    trimmomatic PE -threads 6 -phred33 \
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
fastqc -t 6 "$TRIMMED_DIR"/*.fastq.gz -o "$TRIMMED_QC_DIR"
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
    -p 6 \
    2> "${HOST_ALIGNED_DIR}/${prefix}_bowtie2_report.txt" | \
    samtools sort -@ 6 \
    -o "${HOST_ALIGNED_DIR}/${prefix}_host_aligned_sorted.bam"
done

multiqc "${HOST_ALIGNED_DIR}" -o "${HOST_ALIGNED_DIR}/multiqc_report"

####### 05. Co-Assembly of Metagenomes #######
R1_FILES=""
R2_FILES=""

for r1 in "$NONHOST_DIR"/*_nonhost_1.fastq.gz; do
    r2="${r1/_nonhost_1.fastq.gz/_nonhost_2.fastq.gz}"

    R1_FILES="${R1_FILES:+$R1_FILES,}$r1"
    R2_FILES="${R2_FILES:+$R2_FILES,}$r2"
done

megahit \
    -1 "$R1_FILES" \
    -2 "$R2_FILES" \
    -o "${ASSEMBLY_DIR}/co_assembly" \
    -t 6

####### 06. Quality Control of Coassembly #######
metaquast.py "${ASSEMBLY_DIR}/co_assembly/final.contigs.fa" \
  -o "${ASSEMBLY_DIR}/co_assembly_QC/" \
  -l "MEGAHIT_Metagenomic_Coassembly" \
  --max-ref-num 0 \
  -t 6

####### 07. Discarding Short Contigs #######
megahit_core filterbylen 500 < "${ASSEMBLY_DIR}/co_assembly/final.contigs.fa" > "${ASSEMBLY_DIR}/co_assembly/final.contigs.over500bp.fa"
megahit_core filterbylen 1000 < "${ASSEMBLY_DIR}/co_assembly/final.contigs.fa" > "${ASSEMBLY_DIR}/co_assembly/final.contigs.over1000bp.fa"

####### 08. Alignment to Assembled Metagenome #######
bowtie2-build "${ASSEMBLY_DIR}/co_assembly/final.contigs.over500bp.fa" "${ASSEMBLY_DIR}/co_assembly/final.contigs.over500bp.index"
bowtie2-build "${ASSEMBLY_DIR}/co_assembly/final.contigs.over1000bp.fa" "${ASSEMBLY_DIR}/co_assembly/final.contigs.over1000bp.index"

for r1 in "$NONHOST_DIR"/*_nonhost_1.fastq.gz; do
    r2="${r1/_nonhost_1.fastq.gz/_nonhost_2.fastq.gz}"
    prefix="$(basename "${r1/_nonhost_1.fastq.gz/}")"
    
    echo "Aligning reads to host reference for sample: $prefix"
    
    bowtie2 -x "${ASSEMBLY_DIR}/co_assembly/final.contigs.over500bp.index" \
    -1 "$r1" \
    -2 "$r2" \
    --un-conc-gz "${ALIGNED_DIR}/over500bp/assembly_unaligned/${prefix}_assembly_unaligned_%.fastq.gz" \
    -p 6 \
    2> "${ALIGNED_DIR}/over500bp/${prefix}_bowtie2_report.txt" | \
    samtools sort -@ 6 \
    -o "${ALIGNED_DIR}/over500bp/${prefix}_assembly_aligned_sorted.bam"
done

multiqc "${ALIGNED_DIR}/over500bp" -o "${ALIGNED_DIR}/over500bp/multiqc_report"

####### 09. Assessing Coverage and Binning Contigs #######
jgi_summarize_bam_contig_depths --outputDepth "${ALIGNED_DIR}/over500bp/depth.txt" ${ALIGNED_DIR}/over500bp/*assembly_aligned_sorted.bam

metabat2 -i "${ASSEMBLY_DIR}/co_assembly/final.contigs.over500bp.fa" -a "${ALIGNED_DIR}/depth.txt" -o "${BIN_DIR}/bin"

####### 10. Quality Checking Bins and Mapping Contigs to Bins #######
mkdir "${BIN_DIR}/binning_qc"
checkm2 predict --threads 30 --input "${BIN_DIR}/bin" --output-directory "${BIN_DIR}/binning_qc" --database_path "${PROJECT_DIR}/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd" -x fa

# Get header and determine number of coverage columns
awk 'NR==1 {
    printf "Bin"
    for (i=3; i<=NF; i++)
        printf "\t%s", $i
    printf "\n"
}' "${ALIGNED_DIR}/over500bp/depth.txt" > "${BIN_DIR}/binning_qc/coverage_comparison.tsv"

for bin in "$BIN_DIR"/bin/*.fa; do
    bin_name=$(basename "$bin" .fa)
    # Get contig names in this bin
    grep '^>' "$bin" | sed 's/^>//' | cut -d' ' -f1 > /tmp/bin_contigs.txt
    # Calculate length-weighted coverage for every sample
    awk '
    NR==FNR {
        contig[$1]=1
        next
    }
    FNR==1 {
        next
    }
    ($1 in contig) {
        total_len += $2
        for (i=3; i<=NF; i++)
            cov[i] += $2 * $i
    }
    END {
        printf "%s", "'"$bin_name"'"
        for (i=3; i<=NF; i++) {
            if (total_len > 0)
                printf "\t%.2f", cov[i] / total_len
            else
                printf "\tNA"
        }
        printf "\n"
    }
    ' /tmp/bin_contigs.txt "${ALIGNED_DIR}/over500bp/depth.txt" >> "${BIN_DIR}/binning_qc/coverage_comparison.tsv"
done
rm /tmp/bin_contigs.txt


for f in "${BIN_DIR}"/bin/bin.*.fa; do
    bin=$(basename "$f" .fa)
    grep '^>' "$f" | sed 's/^>//' | cut -d' ' -f1 | \
        awk -v bin="$bin" '{print $1"\t"bin}'
done > "${BIN_DIR}/contig_to_bin.tsv"

####### 11. Taxonomic Annotation with Kraken2 #######
mkdir -p $TAXONOMY_DIR/kraken2_standard_db/trimmed_reads/classified_seqs $TAXONOMY_DIR/kraken2_standard_db/trimmed_reads/unclassified_seqs $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over500bp/classified_seqs $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over500bp/unclassified_seqs $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over1000bp/classified_seqs $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over1000bp/unclassified_seqs $TAXONOMY_DIR/kraken2_standard_db/binned_reads/classified_seqs $TAXONOMY_DIR/kraken2_standard_db/binned_reads/unclassified_seqs

# Annotating trimmed reads
for r1 in "$TRIMMED_DIR"/*_R1_paired.fastq.gz; do
    r2="${r1/_R1_paired.fastq.gz/_R2_paired.fastq.gz}"
    prefix="$(basename "${r1/_R1_paired.fastq.gz/}")"
    
    echo "Kraken2 classification of trimmed reads for sample: $prefix"
    
    kraken2 \
        --db $PROJECT_DIR/kraken2/standard_db \
        --threads 6 \
        --paired \
        --gzip-compressed \
        --report $TAXONOMY_DIR/kraken2_standard_db/trimmed_reads/${prefix}_report.txt \
        --output $TAXONOMY_DIR/kraken2_standard_db/trimmed_reads/${prefix}_output.txt \
        --classified-out $TAXONOMY_DIR/kraken2_standard_db/trimmed_reads/classified_seqs/${prefix}_#.fq \
        --unclassified-out $TAXONOMY_DIR/kraken2_standard_db/trimmed_reads/unclassified_seqs/${prefix}_#.fq \
        ${r1} \
        ${r2} 
done

# Annotating assembled reads
echo "Kraken2 classification of trimmed reads for 500bp co-assembly"
    
kraken2 \
    --db $PROJECT_DIR/kraken2/standard_db \
    --threads 6 \
    --report $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over500bp_report.txt \
    --output $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over500bp_output.txt \
    --classified-out $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over500bp_classified_seqs.fa \
    --unclassified-out $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over500bp_unclassified_seqs.fa \
    $ASSEMBLY_DIR/co_assembly/final.contigs.over500bp.fa

echo "Kraken2 classification of trimmed reads for 1000bp co-assembly"
    
kraken2 \
    --db $PROJECT_DIR/kraken2/standard_db \
    --threads 6 \
    --report $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over1000bp_report.txt \
    --output $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over1000bp_output.txt \
    --classified-out $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over1000bp_classified_seqs.fa \
    --unclassified-out $TAXONOMY_DIR/kraken2_standard_db/assembled_reads/over1000bp_unclassified_seqs.fa \
    $ASSEMBLY_DIR/co_assembly/final.contigs.over1000bp.fa

# Annotating binned sequences
for sample in $BIN_DIR/bin/*.fa; do
    prefix="$(basename "$sample" .fa)"

    echo "Kraken2 classification of binned sequences for bin: $prefix"  

    kraken2 \
    --db $PROJECT_DIR/kraken2/standard_db \
    --threads 6 \
    --report $TAXONOMY_DIR/kraken2_standard_db/binned_reads/${prefix}_report.txt \
    --output $TAXONOMY_DIR/kraken2_standard_db/binned_reads/${prefix}_output.txt \
    --classified-out $TAXONOMY_DIR/kraken2_standard_db/binned_reads/classified_seqs/${prefix}.fa \
    --unclassified-out $TAXONOMY_DIR/kraken2_standard_db/binned_reads/unclassified_seqs/${prefix}.fa \
    ${sample}

done  

# GTDB
mkdir $TAXONOMY_DIR/gtdb-tk_binned_reads
gtdbtk classify_wf \
    --genome_dir $BIN_DIR/bin \
    --out_dir $TAXONOMY_DIR/gtdb-tk_binned_reads \
    --extension fa \
    --cpus 6 \
    --pplacer_cpus 2 \
    --scratch_dir $TAXONOMY_DIR/gtdbtk_scratch

####### 12. Functional Annotation #######

conda activate DRAM
DRAM-setup.py prepare_databases --output_dir DRAM_data --skip-uniref
DRAM.py annotate \
  -i "$BIN_DIR/bin/*.fa" \
  -o DRAM_data/annotation \
  --threads 6 \
  --low_mem_mode


####### 12. Binning Contigs (differently) and Refining Bins #######
