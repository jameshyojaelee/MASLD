#!/bin/bash
#SBATCH --job-name=agp_b1-dbsnp-refresh
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b1_dbsnp_refresh-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b1_dbsnp_refresh-%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT="$1"

echo "host $(hostname)  start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "out  $OUT"
"$PY" "$PROJ/scripts/analysis/alphagenome_program/b1_dbsnp_refresh.py" "$OUT"
echo "end $(date -u +%Y-%m-%dT%H:%M:%SZ)"
