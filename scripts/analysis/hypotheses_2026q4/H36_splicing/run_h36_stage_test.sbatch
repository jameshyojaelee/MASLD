#!/bin/bash
#SBATCH --job-name=H36_stage_test
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H36-20261009T180345Z/logs/slurm-stage_test-%j.out

# H36 step: stage_test. One job, no array. Implements PRESPEC_H36.md.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$REPO/Analysis/Hypotheses/H36-20261009T180345Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H36_splicing
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MASLD_PROJECT_ROOT=$REPO
echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}"
test -f "$OUT/PRESPEC_H36.md"
test -f "$OUT/build/long.rds"
grep -q '^PASS' "$OUT/validation/GATE.txt"   # the real test runs only after the validation gate passed
test ! -e "$OUT/stage_test"
mkdir -p "$OUT/stage_test"
sha256sum "$OUT/PRESPEC_H36.md" "$SC"/h36_dm.R "$SC"/h36_prep.R "$SC"/h36_03_stage_test.R | tee "$OUT/logs/script_sha256_stage_test_${SLURM_JOB_ID:-NA}.txt"
"$RS" "$SC/h36_03_stage_test.R" "$REPO" "$OUT/build" "$OUT/stage_test" 10 0
echo "end $(date -u +%FT%TZ)"
