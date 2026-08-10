#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig5_passport_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig5_passport_%j.err
# ─────────────────────────────────────────────────────────────────────────────
# Fig 5 evidence-passport panels.
#   step 1  (spatial env, pyarrow) reshape the sealed PASS06 candidate -> CSV
#   step 2  (rnaseq  env, ggplot2) render 2 PDFs
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$BASE"
mkdir -p "$BASE/scripts/figures/logs"

echo "=== [1/2] extract passport tables (spatial) ==="
micromamba run -n spatial python -B "$BASE/scripts/figures/fig5_passport_panels.py"

echo "=== [2/2] render panels (rnaseq) ==="
micromamba run -n rnaseq Rscript "$BASE/scripts/figures/fig5_passport_panels.R"

echo "=== done ==="
ls -la "$BASE/figures/main/fig5_convergence/panels/" | grep passport
