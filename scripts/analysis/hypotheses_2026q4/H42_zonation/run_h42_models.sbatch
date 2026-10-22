#!/bin/bash
#SBATCH --job-name=H42_models
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H42-20261009T180511Z/logs/models-%j.out
#SBATCH --error=Analysis/Hypotheses/H42-20261009T180511Z/logs/models-%j.err
# H42 step 2: zonation score versus PNPLA3 I148M dosage, controls and sensitivities (PRESPEC_H42.md plus Addendum 1).
set -euo pipefail
H=Analysis/Hypotheses/H42-20261009T180511Z
GENO=$H/results/geno-20261009T181138Z
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
OUT=$H/results/models-$(date -u +%Y%m%dT%H%M%SZ)
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
"$RS" scripts/analysis/hypotheses_2026q4/H42_zonation/02_zonation_models.R "$GENO" "$OUT"
echo "OUT=$OUT"
