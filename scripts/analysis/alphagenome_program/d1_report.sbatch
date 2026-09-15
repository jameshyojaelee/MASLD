#!/usr/bin/env bash
#SBATCH --job-name=agp_d1_report
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_d1_report_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_d1_report_%j.err
set -euo pipefail
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${OUT:?set OUT}
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
echo "host=$(hostname) start=$(date -u +%FT%TZ)"
"$PY" "$PROJECT/scripts/analysis/alphagenome_program/d1_06_report.py" "$OUT"
"$PY" -m pip freeze > "$OUT/pip_freeze_rnaseq.txt"
echo "end=$(date -u +%FT%TZ)"
