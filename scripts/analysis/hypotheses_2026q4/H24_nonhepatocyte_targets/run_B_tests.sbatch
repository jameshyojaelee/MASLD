#!/bin/bash
#SBATCH --job-name=H24_tests
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# H24 job B: L1/L2 links (05), coloc layers L0/L3/L4/L5 (06), prespecified calls (07).
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_B_tests.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
HERE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H24_nonhepatocyte_targets
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export OMP_NUM_THREADS=4 PATH="$ENV/bin:$PATH"
test -f "$OUT/PRESPEC_H24.md"
N_PERM=$($ENV/bin/python -I -c "import json,sys; print(json.load(open(sys.argv[1]))['n_perm'])" "$OUT/inventory_v2/inventory.json")
echo "N_PERM=$N_PERM"
$ENV/bin/Rscript "$HERE/05_links.R" "$OUT" "$N_PERM"
$ENV/bin/Rscript "$HERE/06_qtl_coloc.R" "$OUT"
$ENV/bin/Rscript "$HERE/07_calls.R" "$OUT"
echo "H24 job B done"
