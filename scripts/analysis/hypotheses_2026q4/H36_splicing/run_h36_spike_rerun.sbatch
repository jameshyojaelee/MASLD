#!/bin/bash
#SBATCH --job-name=H36_spike_rerun
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H36-20261009T180345Z/logs/slurm-spike_rerun-%j.out

# H36: re-run of validation part (iii) with the corrected true dPSI_0to4. One job, no array.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$REPO/Analysis/Hypotheses/H36-20261009T180345Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H36_splicing
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}"
test ! -e "$OUT/validation_spike_rerun"
mkdir -p "$OUT/validation_spike_rerun"
sha256sum "$SC"/h36_02b_spike_rerun.R | tee "$OUT/logs/script_sha256_spike_rerun_${SLURM_JOB_ID:-NA}.txt"
"$RS" "$SC/h36_02b_spike_rerun.R" "$REPO" "$OUT/build" "$OUT/validation_spike_rerun"
echo "end $(date -u +%FT%TZ)"
