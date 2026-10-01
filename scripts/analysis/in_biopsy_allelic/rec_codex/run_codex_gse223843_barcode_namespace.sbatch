#!/usr/bin/env bash
#SBATCH --job-name=1592
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1592_barcode_namespace_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1592_barcode_namespace_%j.err
set -euo pipefail
umask 077
ulimit -c 0
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
exec /usr/bin/python3 "$TASK_ROOT/scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse223843_barcode_namespace.py" \
    --expected-script-sha256 cb269cd1f6e833bec756603e6d3bf60655cd5b140f8add4b0b3c80b5d252745c
