#!/bin/bash
#SBATCH --job-name=H23_msa_fm
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --output=Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.out
#SBATCH --error=Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.err
# H23: EUR-only SuSiE-RSS refits at MESuSiE loci, then endpoints and tests.
# Usage: sbatch [--cpus-per-task=N --mem=M --job-name=H23_smoke] run_h23.sbatch <smoke|full>
set -eo pipefail
MODE="${1:?mode smoke|full}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
H=$ROOT/Analysis/Hypotheses/H23-20261009T140416Z
S=$ROOT/scripts/analysis/hypotheses_2026q4/H23_multiancestry_finemapping
cd "$ROOT"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
export H23_CORES="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # one thread per forked worker
[ "$MODE" = smoke ] && export H23_MIN_UNITS=2 || true   # smoke only; the prespecified minimum is 10
TS=$(date -u +%Y%m%dT%H%M%SZ)
OUT="$H/results/${MODE}-${TS}"
Rscript "$S/h23_fit_units.R" "$OUT/fits" "$([ "$MODE" = smoke ] && echo smoke || echo full)"
Rscript "$S/h23_analyze.R" "$OUT/fits" "$OUT/analysis"
pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
echo "done $OUT"
