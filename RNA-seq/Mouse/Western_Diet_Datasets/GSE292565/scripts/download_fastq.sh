#!/bin/bash
#SBATCH --job-name=sra_GSE292565
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=../logs/download_%A_%a.out
#SBATCH --error=../logs/download_%A_%a.err

# Download FASTQs for GSE292565 (FFC diet, paired-end DNBSEQ-T7)
# Usage: sbatch --array=1-24 download_fastq.sh

set -euo pipefail

module load SRA-Toolkit/3.2.0-gompi-2023b
module load pigz/2.7-GCCcore-12.2.0

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE292565"
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BASE}/metadata/srr_list.txt")

if [[ -z "${SRR}" ]]; then
    echo "ERROR: No SRR for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "[$(date)] Downloading ${SRR} (task ${SLURM_ARRAY_TASK_ID})"

cd "${BASE}/fastq"

# Prefetch to local cache then extract
prefetch --force all "${SRR}" -O .
fasterq-dump --split-files --threads 4 "${SRR}/${SRR}.sra" -O .
rm -rf "${SRR}"  # clean up .sra cache

# Compress with pigz
pigz -p 4 "${SRR}"*.fastq

echo "[$(date)] Done: ${SRR}"
ls -lh "${SRR}"*.fastq.gz
