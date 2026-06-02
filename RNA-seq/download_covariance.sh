#!/bin/bash
#SBATCH --job-name=dl_cov
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/dl_cov_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/dl_cov_%j.err

GWAS_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data"
TMP_TAR="/tmp/mashr_eqtl_${SLURM_JOB_ID}.tar"

echo "Starting download at $(date)"
# Use wget with retries
wget -t 3 -c -O "$TMP_TAR" "https://zenodo.org/records/3518299/files/mashr_eqtl.tar?download=1"

if [ $? -eq 0 ]; then
    echo "Download successful. Starting extraction."
    tar -xf "$TMP_TAR" -C /tmp/ eqtl/mashr/mashr_Liver.txt.gz
    
    if [ -f "/tmp/eqtl/mashr/mashr_Liver.txt.gz" ]; then
        echo "Extraction successful. Moving to $GWAS_DIR"
        mv "/tmp/eqtl/mashr/mashr_Liver.txt.gz" "$GWAS_DIR/mashr_Liver.txt.gz"
        echo "Done."
    else
        echo "Error: extracted file not found in expected location."
        exit 1
    fi
else
    echo "Error: Download failed."
    exit 1
fi

# Cleanup
rm -rf /tmp/eqtl "$TMP_TAR"
echo "Finished at $(date)"
