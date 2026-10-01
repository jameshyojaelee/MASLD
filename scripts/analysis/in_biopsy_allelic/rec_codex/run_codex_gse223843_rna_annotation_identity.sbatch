#!/usr/bin/env bash
#SBATCH --job-name=1592
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1592_rna_annotation_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1592_rna_annotation_%j.err
set -euo pipefail
umask 077
ulimit -c 0
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
exec /usr/bin/python3 "$TASK_ROOT/scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse223843_rna_annotation_identity.py" \
    --expected-script-sha256 8620b8365044673de0c39ed2e1fe9875b88057c93860c57b4fb04ec484a8a804
