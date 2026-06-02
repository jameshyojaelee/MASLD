#!/bin/bash
#SBATCH --job-name=347c_sh_dissoc
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/347c_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/347c_%j.err

set -euo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE"
mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs

export MASLD_PROJECT_ROOT="$BASE"

echo "[$(date)] 347c: SH-vs-Healthy dissociation sensitivity"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/347c_gse136103_holdout_sensitivity_sh.R

echo "[$(date)] 349c: re-annotate dedup top-30 against combined dissoc table"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/349c_annotate_dissoc_status.R

echo "[$(date)] Done"
