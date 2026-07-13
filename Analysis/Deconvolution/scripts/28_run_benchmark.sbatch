#!/bin/bash
#SBATCH --job-name=rectangle
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/scripts/logs/28_benchmark_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/scripts/logs/28_benchmark_%j.log
#
# Pseudobulk ground-truth deconvolution benchmark: Rectangle vs MuSiC vs InstaPrism on
# simulated bulk with KNOWN cell-type proportions. Chains 28a -> 28b -> 28c -> 28d, each
# in its native env (method<->modality split). Fallback if cpu is busy:
#   sbatch --partition=bigmem --qos=interactive --mem=500G 28_run_benchmark.sbatch
#   sbatch --partition=io     --qos=interactive --mem=96G  28_run_benchmark.sbatch
set -euo pipefail

ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPTS="$ROOT/Analysis/Deconvolution/scripts"
export BENCH_OUT_DIR="${BENCH_OUT_DIR:-$ROOT/Analysis/Deconvolution/rectangle_comparison/pseudobulk}"
export BENCH_MIN_CELLS="${BENCH_MIN_CELLS:-200}"
export BENCH_N_FOLDS="${BENCH_N_FOLDS:-5}"
export BENCH_MAX_DONORS="${BENCH_MAX_DONORS:-0}"
export BENCH_RECT_OPTIMIZE="${BENCH_RECT_OPTIMIZE:-1}"

PY_RECT="$ROOT/RNA-seq/.mamba/rectangle/bin/python"
R_MUSIC="$ROOT/RNA-seq/.mamba/music_deconv/bin/R"   # bin/R (the env's Rscript shim is broken)

echo "=== benchmark config ==="
echo "OUT_DIR=$BENCH_OUT_DIR  MIN_CELLS=$BENCH_MIN_CELLS  N_FOLDS=$BENCH_N_FOLDS"
echo "MAX_DONORS=$BENCH_MAX_DONORS  RECT_OPTIMIZE=$BENCH_RECT_OPTIMIZE  CPUS=${SLURM_CPUS_PER_TASK:-NA}"
mkdir -p "$BENCH_OUT_DIR"

echo "=== [1/4] 28a pseudobulk prepare (rectangle env) ==="
srun --ntasks=1 "$PY_RECT" "$SCRIPTS/28a_pseudobulk_prepare.py"

echo "=== [2/4] 28b Rectangle (rectangle env) ==="
srun --ntasks=1 "$PY_RECT" "$SCRIPTS/28b_run_rectangle_benchmark.py"

echo "=== [3/4] 28c MuSiC + InstaPrism (music_deconv env) ==="
srun --ntasks=1 "$R_MUSIC" --vanilla --no-echo -f "$SCRIPTS/28c_run_r_benchmark.R"

echo "=== [4/4] 28d score (rectangle env) ==="
srun --ntasks=1 "$PY_RECT" "$SCRIPTS/28d_score_benchmark.py"

echo "=== benchmark complete ==="
