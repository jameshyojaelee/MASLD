#!/bin/bash
#SBATCH --job-name=agp_e1-bridge
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e1_bridge_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e1_bridge_%j.err
set -euo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=$PROJ/GWAS/finemapping/results/alphagenome_program/e1-bridge-20260914T194145Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1

mkdir -p "$OUT"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} out=$OUT"
"$PY" "$PROJ/scripts/analysis/alphagenome_program/e1_bridge.py" --outdir "$OUT" --threads 8
cp "$PROJ/scripts/analysis/alphagenome_program/e1_bridge.py" "$OUT/e1_bridge.py.snapshot"
cp "$PROJ/scripts/analysis/alphagenome_program/e1_bridge.sbatch" "$OUT/e1_bridge.sbatch.snapshot"
echo "COMPLETE"
