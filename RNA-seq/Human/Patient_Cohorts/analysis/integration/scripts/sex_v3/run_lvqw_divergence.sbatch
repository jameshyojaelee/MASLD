#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x_lvqwdiv_%j.out
#SBATCH --error=%x_lvqwdiv_%j.out
# ---------------------------------------------------------------------------
# LVQW divergence report — runs AFTER 06 (chained afterok). Compares staging
# LVQW sex calls vs canonical dream. Survives harness restarts (SLURM-native).
# NOTE: NO `set -u` (conda activate.d ADDR2LINE unbound-var). -e + pipefail only.
# ---------------------------------------------------------------------------
set -eo pipefail
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

cd "$SCRIPTS"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== LVQW divergence report ==="
echo "host: $(hostname)  started: $(date)"
Rscript lvqw_divergence_report.R
echo "finished divergence: $(date)"
