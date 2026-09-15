#!/bin/bash
#SBATCH --job-name=agp_b1-response-layer-aux
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b1-response-layer-aux_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b1-response-layer-aux_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
S=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/b1-response-layer-20260914T195729Z
$PY $S/b1_mpra_rebuild.py "$OUT"
$PY $S/b1_direction_track_panel.py "$OUT"
$PY $S/b1_checks.py "$OUT"
