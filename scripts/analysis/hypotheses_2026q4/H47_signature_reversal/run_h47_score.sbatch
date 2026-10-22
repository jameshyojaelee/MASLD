#!/bin/bash
#SBATCH --job-name=H47score
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# H47 step B: unit checks, then reversal scoring, class tests, null, sensitivity.
# Usage: sbatch --output=<H47>/logs/score-%j.log run_h47_score.sbatch <H47_dir> <extract_dir> <run_dir>
# Dry run on synthetic data: H47_DRYRUN=1 (builds a synthetic extract in <extract_dir>, 20 null draws).
set -euo pipefail
H47="${1:?H47 dir required}"
EXT="${2:?extract dir required}"
OUT="${3:?run dir required}"
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H47_signature_reversal
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
DEG=$ROOT/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv
META=$ROOT/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 R_LIBS_USER= PYTHONHASHSEED=0
test -s "$H47/PRESPEC_H47.md"
test ! -e "$OUT"
cd "$SCRIPTS"
"$ENV/bin/python3" -I test_h47_score.py
if [ "${H47_DRYRUN:-0}" = "1" ]; then
  export H47_N_NULL=20
  "$ENV/bin/python3" -I dryrun_synthetic_extract.py "$H47/downloads_untrusted" "$EXT"
  "$ENV/bin/Rscript" --vanilla h47_map_landmarks.R "$H47/downloads_untrusted/GSE92742_Broad_LINCS_gene_info.txt.gz" "$DEG" "$META" "$EXT"
fi
test -s "$EXT/landmark_map.tsv"
mkdir -p "$(dirname "$OUT")"
"$ENV/bin/python3" -I h47_score.py "$H47" "$EXT" "$DEG" "$OUT"
sha256sum h47_score.py test_h47_score.py h47_extract.py h47_map_landmarks.R run_h47_score.sbatch run_h47_fetch.sbatch > "$OUT/script_sha256.txt"
sha256sum "$DEG" "$EXT"/*.npy "$EXT"/*.tsv "$H47/PRESPEC_H47.md" > "$OUT/input_sha256.txt"
cp h47_score.py "$OUT/h47_score.as_run.py"
"$ENV/bin/python3" -I -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
"$ENV/bin/python3" --version > "$OUT/python_version.txt"
cp "$EXT/sessionInfo_map.txt" "$OUT/" 2>/dev/null || true
echo "[$(date -u +%FT%TZ)] done"
