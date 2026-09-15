#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_bench_eval
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z/logs/bench_eval-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z
export PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
echo "== caqtl (kill criterion: Atlas ATAC must reproduce the archived API auROC)"
$PY 33_benchmark_eval.py caqtl
echo "== sqtl"
$PY 33_benchmark_eval.py sqtl
echo BENCH_EVAL_DONE
