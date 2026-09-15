#!/bin/bash
#SBATCH --job-name=agp_b1-response-layer
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b1-response-layer_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b1-response-layer_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/b1-response-layer-20260914T195729Z
echo "host $(hostname) start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "out $OUT"
$PY /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program/b1_response_layer.py "$OUT"
echo "done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
