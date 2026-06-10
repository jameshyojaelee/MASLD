#!/bin/bash
#SBATCH --job-name=Rexact
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/degx_full_battery_rerun_%j.log
# Fresh re-run of the FULL degx R-exact method battery (all ~19 DE variants) on the
# identical 5-cohort disease_vs_control data, into the project, to compare against
# the stored canonical degx run (~/degx/runs/Rexact/).
# cpu is default; when congested override to io/bigmem fast lane:
#   sbatch --partition=io --qos=interactive --mem=96G run_degx_full_battery_rerun.sh

set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
DEGX_HOME=/gpfs/commons/home/jameslee/degx
OUT="$PROJECT_ROOT/figures/supplementary/figS_methods_validation/disease_signature_sweep/deg_tables/full_battery"

mkdir -p "$PROJECT_ROOT/logs" "$OUT"
cd "$DEGX_HOME"   # engine sources R/corrections.R relative to getwd()

echo "[$(date)] degx full-battery re-run on $(hostname)  JOB=${SLURM_JOB_ID:-local}  CPUs=${SLURM_CPUS_PER_TASK:-8}"
echo "out_dir = $OUT"

micromamba run -n rnaseq Rscript R/run_methods_real.R \
  --contrast disease_vs_control \
  --data_dir "$DEGX_HOME/data/masld" \
  --out_dir  "$OUT"

echo "[$(date)] Done. Tables in $OUT"
