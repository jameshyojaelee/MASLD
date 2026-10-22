#!/bin/bash
#SBATCH --job-name=H36_integrate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H36-20261009T180345Z/logs/slurm-integrate-%j.out

# H36 step: calls, calibration, consistency, S1, leave-one-cohort-out, concordance with sQTL colocalization. One job, no array.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$REPO/Analysis/Hypotheses/H36-20261009T180345Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H36_splicing
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}"
test -f "$OUT/stage_test_v2/stage_test_observed.tsv.gz"
test -f "$OUT/sqtl_coloc/sqtl_coloc_locus_intron.tsv.gz"
test ! -e "$OUT/results/integrate_v2"
mkdir -p "$OUT/results/integrate_v2"
sha256sum "$OUT/PRESPEC_H36.md" "$SC"/h36_dm.R "$SC"/h36_prep.R "$SC"/h36_05_integrate.R | tee "$OUT/logs/script_sha256_integrate_${SLURM_JOB_ID:-NA}.txt"
"$RS" "$SC/h36_05_integrate.R" "$REPO" "$OUT" "$OUT/results/integrate_v2" stage_test_v2 0
echo "end $(date -u +%FT%TZ)"
