#!/bin/bash
#SBATCH --job-name=download_hra007511
#SBATCH --output=scripts/data_retrieval/logs/%x_%A.out
#SBATCH --error=scripts/data_retrieval/logs/%x_%A.err
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00

echo "Starting CNCB NGDC Fastq FTP Download for HRA007511..."
mkdir -p data/Spatial/HRA007511/fastq

# Download the dataset cleanly parsing directory structures
wget -c -r --no-parent -nH --cut-dirs=2 "ftp://download.big.ac.cn/gsa-human/HRA007511/" -P data/Spatial/HRA007511/fastq/

if [ $? -eq 0 ]; then
    echo "SUCCESS: HRA007511 fastqs downloaded successfully."
else
    echo "ERROR: wget encountered a failure pulling HRA007511."
    exit 1
fi
