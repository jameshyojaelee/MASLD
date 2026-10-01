#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
# Model A-final B-MODEL jobs (synthetic inputs only; no MASLD allele data).
# Usage: sbatch --job-name=<12xx> -o <out>/logs/%j.out run_model.sbatch <out> <mode> [args]
#   mode pytest  -> pytest on tests/, log in <out>/logs
#   mode fixture -> run_fixture.py on the B-FIX fixture into <out>
#   mode py      -> any python script in this directory with args
set -euo pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
HERE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/model
FIX=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z
OUT=$1
MODE=$2
shift 2
cd "$HERE"
case "$MODE" in
  pytest) $PY -m pytest -q -p no:cacheprovider tests/ "$@" ;;
  fixture)
    mkdir -p "$OUT"
    sha256sum "$HERE"/*.py "$HERE"/tests/*.py "$HERE"/run_model.sbatch > "$OUT/sha256_code.txt"
    $PY -m pip freeze > "$OUT/pip_freeze.txt"
    $PY run_fixture.py --fixture "$FIX" --out "$OUT" "$@" ;;
  py) $PY "$@" ;;
esac
