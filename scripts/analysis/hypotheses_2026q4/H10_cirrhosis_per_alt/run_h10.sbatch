#!/bin/bash
#SBATCH --job-name=H10_cirr_per_alt
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H10-20261007T210859Z/slurm-%j.out

set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$REPO/Analysis/Hypotheses/H10-20261007T210859Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
test -f "$OUT/PRESPEC_H10.md"
sha256sum "$OUT/PRESPEC_H10.md" "$REPO/scripts/analysis/hypotheses_2026q4/H10_cirrhosis_per_alt/h10_cirrhosis_per_alt.py"
"$PY" "$REPO/scripts/analysis/hypotheses_2026q4/H10_cirrhosis_per_alt/h10_cirrhosis_per_alt.py" "$REPO" "$OUT"
