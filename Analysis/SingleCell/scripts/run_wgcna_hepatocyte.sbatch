#!/bin/bash
#SBATCH --job-name=wgcna_hep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/wgcna_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/wgcna_%j.err

set -euo pipefail

SC_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell"
mkdir -p "${SC_DIR}/logs"

echo "=== WGCNA Hepatocyte Analysis ==="
echo "Start: $(date)"

module purge 2>/dev/null || true
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript "${SC_DIR}/scripts/fig2_wgcna_hepatocyte.R"

echo "=== Done: $(date) ==="
ls -la "${SC_DIR}/results_gpu_v2/fig2_data/wgcna_"*.csv
