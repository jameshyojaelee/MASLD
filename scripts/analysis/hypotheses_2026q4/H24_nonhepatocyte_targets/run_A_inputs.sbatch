#!/bin/bash
#SBATCH --job-name=H24_inputs
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# H24 job A: loci, pseudobulk inputs and test-universe inventory. No outcome.
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_A_inputs.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
HERE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H24_nonhepatocyte_targets
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
$ENV/bin/python -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
$ENV/bin/python -I "$HERE/01_define_loci.py" "$OUT"
$ENV/bin/python -I "$HERE/02_prepare_gse296875.py" "$OUT"
$ENV/bin/python -I "$HERE/03_prepare_gse244832.py" "$OUT"
$ENV/bin/Rscript "$HERE/04_inventory.R" "$OUT"
echo "H24 job A done"
