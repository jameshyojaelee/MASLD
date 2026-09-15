#!/bin/bash
#SBATCH --job-name=agp_e3-reporter-context-probe
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e3_probe_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e3_probe_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT="$1"
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python

echo "host $(hostname) start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
"$PY" "$PROJ/scripts/analysis/alphagenome_program/e3_score_api.py" "$OUT" --probe --max-calls 8
echo "end $(date -u +%Y-%m-%dT%H:%M:%SZ)"
