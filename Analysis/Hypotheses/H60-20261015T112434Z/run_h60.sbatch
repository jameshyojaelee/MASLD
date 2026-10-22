#!/bin/bash
#SBATCH --job-name=6000
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z/logs/slurm-%x-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
JOB=$ROOT/Analysis/Hypotheses/H60-20261015T112434Z
TAG=${H60_TAG:-final}
OUT=$JOB/out-$TAG-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "$OUT"; mkdir -p "$OUT/tmp"
export TMPDIR=$OUT/tmp XDG_CACHE_HOME=$OUT/tmp MPLCONFIGDIR=$OUT/tmp
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
cd "$ROOT"
echo "OUT=$OUT host=$(hostname) start=$(date -u)"
# 1 input hashes (stops on mismatch)
$ENV/python $JOB/h60_verify_inputs.py "$OUT"
# 2 R0 reproduction and R1 hepatocyte-adjusted refit
$ENV/Rscript $JOB/h60_hepadj_refit.R "$OUT"
# 3 diagnostics, controls, stop rules, then primary, qualifier, fixed-sequence step, descriptives
$ENV/python $JOB/h60_signed_twas_de.py "$OUT"
# 4 environment record
$ENV/python -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
echo "done $(date -u)"
