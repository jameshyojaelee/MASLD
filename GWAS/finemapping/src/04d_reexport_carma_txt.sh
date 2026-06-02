#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=reexport_carma
#SBATCH --output=logs/reexport_carma_%j.out
#SBATCH --error=logs/reexport_carma_%j.err

# Re-export CARMA .txt.gz files from .rds backups
# Fixes 0-byte files caused by data.table::fwrite gzip compression bug
# Usage: sbatch src/04d_reexport_carma_txt.sh

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping

echo "Starting CARMA re-export: $(date)"
Rscript src/04d_reexport_carma_txt.R
echo "Finished CARMA re-export: $(date)"
