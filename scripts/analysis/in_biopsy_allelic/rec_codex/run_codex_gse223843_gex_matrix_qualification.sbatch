#!/usr/bin/env bash
#SBATCH --job-name=1592
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1592_gex_mtx_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1592_gex_mtx_%j.err
set -euo pipefail
umask 077
ulimit -c 0
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
exec /usr/bin/python3 "$TASK_ROOT/scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse223843_gex_matrix_qualification.py" \
    --expected-script-sha256 e3482519319ec45db624a2bb1d92eda366a2cda122e2af7d145f0c67e856186a
