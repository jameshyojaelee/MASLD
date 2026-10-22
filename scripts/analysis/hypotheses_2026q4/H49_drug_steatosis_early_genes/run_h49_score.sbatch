#!/bin/bash
#SBATCH --job-name=H49score
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# H49 step C: unit checks, then scoring.
# Usage: sbatch --output=<H49>/logs/score-%j.log run_h49_score.sbatch <H49_dir> <per_compound_dir> <run_dir>
# Dry run on synthetic inputs: H49_DRYRUN=1 (per_compound_dir is ignored; 1,000 null draws).
set -euo pipefail
H49="${1:?H49 dir required}"
PCD="${2:?per-compound dir required}"
OUT="${3:?run dir required}"
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H49_drug_steatosis_early_genes
H47X=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H47-20261008T214453Z/extract-20261008T215141Z
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 R_LIBS_USER= PYTHONHASHSEED=0
test -s "$H49/PRESPEC_H49.md"
test ! -e "$OUT"
cd "$SCRIPTS"
"$ENV/bin/Rscript" --vanilla test_h49_lib.R
ANN="$H49/downloads_untrusted_annotation/GPL570.annot.gz"
S5="$H49/downloads_untrusted_zhao2026/elife-109534-supp5-v1.xlsx"
S9="$H49/downloads_untrusted_zhao2026/elife-109534-supp9-v1.xlsx"
DILI="$H49/downloads_untrusted_dilirank/DILIrank-List.xlsx"
if [ "${H49_DRYRUN:-0}" = "1" ]; then
  export H49_N_NULL=1000
  SYN="$OUT-inputs"
  "$ENV/bin/Rscript" --vanilla dryrun_synthetic.R "$ANN" "$S5" "$H47X/landmark_order.tsv" "$SYN"
  PCD="$SYN/per_compound"
  LINCS="$SYN/lincs_phh_synthetic.tsv.gz"
else
  mkdir -p "$OUT-inputs"
  LINCS="$OUT-inputs/lincs_phh_z.tsv.gz"
  "$ENV/bin/python3" -I h49_lincs_export.py "$H47X" "$LINCS"
fi
"$ENV/bin/Rscript" --vanilla h49_score.R "$PCD" "$ANN" "$S5" "$S9" "$LINCS" "$DILI" "$OUT"
sha256sum h49_score.R h49_lib.R h49_lincs_export.py test_h49_lib.R > "$OUT/script_sha256.txt"
"$ENV/bin/python3" -I -m pip freeze > "$OUT/pip_freeze.txt"
