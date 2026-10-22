#!/bin/bash
#SBATCH --job-name=H41
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
# H41: extract HEPG2 shRNA landmark signatures (GSE92742 level 5), map landmarks,
# score knockdown consensus signatures against the refit disease signature.
# Usage: sbatch --output=<H41>/logs/h41-%j.log run_h41.sbatch <H41_dir> <extract_dir> <run_dir>
# Dry run on synthetic z-scores: H41_DRYRUN=1 (H41_N_NULL may lower the draw count).
set -euo pipefail
H41="${1:?H41 dir required}"
EXT="${2:?extract dir required}"
OUT="${3:?run dir required}"
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
SCRIPTS=$ROOT/scripts/analysis/hypotheses_2026q4/H41_identity_tf_knockdown
MAPR=$ROOT/scripts/analysis/hypotheses_2026q4/H47_signature_reversal/h47_map_landmarks.R
H47=$ROOT/Analysis/Hypotheses/H47-20261008T214453Z
DL=$H47/downloads_untrusted
DEC=$H47/downloads_untrusted_decompressed
GCTX=$DEC/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
GCTX_SHA=b293f3fb7c2298a60526de727e5400d8400af4b77a23c4ed2116f86199fb45e8
DEG=$ROOT/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv
META=$ROOT/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 R_LIBS_USER= PYTHONHASHSEED=0
test -s "$H41/PRESPEC_H41.md"
test -s "$H41/negative_control_genes.tsv"
test ! -e "$EXT"
test ! -e "$OUT"
cd "$SCRIPTS"
echo "[$(date -u +%FT%TZ)] unit tests"
"$ENV/bin/python3" -I test_h41_score.py
if [ "${H41_DRYRUN:-0}" = "1" ]; then
  "$ENV/bin/python3" -I h41_dryrun_synthetic_extract.py "$DL" "$EXT"
else
  echo "[$(date -u +%FT%TZ)] S0 sha256 check of the level 5 GCTX"
  echo "$GCTX_SHA  $GCTX" | sha256sum -c -
  echo "[$(date -u +%FT%TZ)] extracting"
  "$ENV/bin/python3" -I h41_extract.py "$DL" "$DEC" "$EXT"
fi
echo "[$(date -u +%FT%TZ)] mapping landmarks"
"$ENV/bin/Rscript" --vanilla "$MAPR" "$DL/GSE92742_Broad_LINCS_gene_info.txt.gz" "$DEG" "$META" "$EXT"
if [ "${H41_DRYRUN:-0}" = "1" ] && [ "${H41_DRYRUN_PLANT:-0}" = "1" ]; then
  "$ENV/bin/python3" -I h41_dryrun_plant.py "$EXT" "$DEG"
fi
echo "[$(date -u +%FT%TZ)] scoring"
mkdir -p "$(dirname "$OUT")"
"$ENV/bin/python3" -I h41_score.py "$H41" "$EXT" "$DEG" "$OUT"
sha256sum h41_score.py h41_extract.py h41_draw_negative_controls.py h41_dryrun_synthetic_extract.py \
  h41_dryrun_plant.py test_h41_score.py run_h41.sbatch "$MAPR" > "$OUT/script_sha256.txt"
sha256sum "$DEG" "$H41/PRESPEC_H41.md" "$H41/negative_control_genes.tsv" "$EXT"/*.npy "$EXT"/*.tsv \
  "$DL/GSE92742_Broad_LINCS_sig_info.txt.gz" "$DL/GSE92742_Broad_LINCS_gene_info.txt.gz" > "$OUT/input_sha256.txt"
cp h41_score.py "$OUT/h41_score.as_run.py"
"$ENV/bin/python3" -I -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
"$ENV/bin/python3" --version > "$OUT/python_version.txt"
cp "$EXT/sessionInfo_map.txt" "$OUT/sessionInfo_R_map.txt"
echo "[$(date -u +%FT%TZ)] done"
