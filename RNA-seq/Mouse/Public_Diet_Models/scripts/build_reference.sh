#!/bin/bash
#SBATCH --job-name=build_vm38_ref
#SBATCH --output=logs/ref_build_%j.out
#SBATCH --error=logs/ref_build_%j.err
#SBATCH --partition=bigmem
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00

# Setup directories
ROOT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
REF_RAW="${ROOT_DIR}/reference/raw"
REF_INDEX="${ROOT_DIR}/reference/star_vm38"
LOGS="${ROOT_DIR}/logs"

mkdir -p "$REF_RAW" "$REF_INDEX" "$LOGS"

# URLs
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz"
FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz"

# Activation (Ensure STAR is available)
source ~/.bashrc
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl

# 1. Download Files
echo "Downloading Reference Files..."
cd "$REF_RAW"

if [ ! -f "gencode.vM38.annotation.gtf" ]; then
    wget -O gencode.vM38.annotation.gtf.gz "$GTF_URL"
    gunzip gencode.vM38.annotation.gtf.gz
fi

if [ ! -f "GRCm39.primary_assembly.genome.fa" ]; then
    wget -O GRCm39.primary_assembly.genome.fa.gz "$FASTA_URL"
    gunzip GRCm39.primary_assembly.genome.fa.gz
fi

# 2. Build STAR Index
echo "Building STAR Index..."
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir "$REF_INDEX" \
     --genomeFastaFiles "$REF_RAW/GRCm39.primary_assembly.genome.fa" \
     --sjdbGTFfile "$REF_RAW/gencode.vM38.annotation.gtf" \
     --sjdbOverhang 149 \
     --limitGenomeGenerateRAM 60000000000

echo "Done."
