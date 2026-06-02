#!/bin/bash
#SBATCH --job-name=download_pxd052937
#SBATCH --output=scripts/data_retrieval/logs/%x_%A.out
#SBATCH --error=scripts/data_retrieval/logs/%x_%A.err
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00

echo "Starting PXD052937 Proteomics FTP Download..."
mkdir -p data/Proteomics/PXD052937

# Download the search (processed) matrices recursively
wget -c -r --no-parent -nH --cut-dirs=3 "ftp://massive-ftp.ucsd.edu/v08/MSV000094959/search/" -P data/Proteomics/PXD052937/

if [ $? -eq 0 ]; then
    echo "SUCCESS: PXD052937 matrices downloaded successfully."
else
    echo "ERROR: wget encountered a failure."
    exit 1
fi
