#!/usr/bin/env bash
#SBATCH --job-name=1594
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1594_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1594_%j.err
set -euo pipefail
umask 027
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
exec /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$TASK_ROOT/scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse296875_compare_native_common_counts.py" \
  --expected-helper-sha 81545a7690cf1e66d2358f99ae55da9ff2605496545e20826f4726225abf6af8
