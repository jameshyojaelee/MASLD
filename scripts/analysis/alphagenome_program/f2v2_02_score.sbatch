#!/usr/bin/env bash
#SBATCH --job-name=agp_f2-haplotype-v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f2-haplotype-v2-%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export F2V2_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/f2-haplotype-v2-20260915T063000Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program
"$PY" -m pip freeze > "$F2V2_OUT_ROOT/pip_freeze.txt" 2>/dev/null || true
"$PY" f2v2_01_draw_genespan_null.py
"$PY" f2v2_02_score_genespan_null.py
echo F2V2_SCORE_DONE
