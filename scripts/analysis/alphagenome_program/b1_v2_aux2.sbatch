#!/bin/bash
#SBATCH --job-name=agp_b1-response-layer-v2-dict
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b1-response-layer-v2-dict_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b1-response-layer-v2-dict_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
S=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program
R=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program
OUT=$R/b1-response-layer-v2-20260914T204700Z
echo "host $(hostname) start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
$PY $S/b1_v2_column_dictionary.py "$OUT"
$PY $S/b1_v2_compare_to_wave1.py "$R/b1-response-layer-20260914T195729Z" "$OUT" "$OUT/numbers_moved_vs_wave1.tsv"
echo "done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
