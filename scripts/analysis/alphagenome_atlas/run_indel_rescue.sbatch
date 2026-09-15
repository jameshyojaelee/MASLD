#!/usr/bin/env bash
#SBATCH --job-name=ag_indel_rescue
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ag_indel_rescue_%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT="$1"
export PYTHONNOUSERSITE=1
export AGA_INDEL_MIN_MASS=${AGA_INDEL_MIN_MASS:-0.05}
export AGA_INDEL_CONTROLS=${AGA_INDEL_CONTROLS:-4}
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd "$MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas"
$PY 63_indel_rescue.py
$PY 64_indel_vs_snv_paired.py
echo "INDEL_RESCUE_DONE"
