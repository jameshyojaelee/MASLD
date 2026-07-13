#!/bin/bash -l
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu,io
#SBATCH --job-name=aggmvp
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/aggmvp_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/aggmvp_%j.err
set -eo pipefail
FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq; set -u
echo "[aggmvp] SuSiEx aggregation $(date)"
Rscript src/11_aggregate_susiex_mvp.R
echo "[aggmvp] meSuSiE aggregation $(date)"
Rscript src/11b_aggregate_mesusie_mvp.R
echo "[aggmvp] done $(date)"
