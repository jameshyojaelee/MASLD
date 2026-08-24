#!/usr/bin/env bash
set -euo pipefail

: "${MASLD_PROJECT_ROOT:?MASLD_PROJECT_ROOT is required}"
: "${FIG2_SUPP_OUT_DIR:?FIG2_SUPP_OUT_DIR is required}"

mkdir -p "$FIG2_SUPP_OUT_DIR"
cd "$MASLD_PROJECT_ROOT"
printf 'job_id\tstarted_utc\n%s\t%s\n' \
  "${SLURM_JOB_ID:-not_slurm}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  > "$FIG2_SUPP_OUT_DIR/PREFLIGHT_STARTED.tsv"

canonical_dir="figures/main/fig2_genetics/panels"
find "$canonical_dir" -maxdepth 1 -type f -regextype posix-extended \
  -regex '.*/FigS2[J-R].*' -print0 | sort -z | xargs -0 -r sha256sum \
  > "$FIG2_SUPP_OUT_DIR/canonical_before.sha256"

RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

"$RS" scripts/figures/fig2_crossancestry_concentration.R \
  > "$FIG2_SUPP_OUT_DIR/FigS2J.log" 2>&1
if "$RS" scripts/figures/figS_cs_regulatory.R \
  > "$FIG2_SUPP_OUT_DIR/FigS2K.log" 2>&1; then
  exit 41
fi
grep -q 'FigS2K retired' "$FIG2_SUPP_OUT_DIR/FigS2K.log"
"$RS" scripts/figures/fig2_f2rl1_spotlight.R \
  > "$FIG2_SUPP_OUT_DIR/FigS2L.log" 2>&1
"$RS" scripts/figures/figS_mesusie_shared.R \
  > "$FIG2_SUPP_OUT_DIR/FigS2M.log" 2>&1
if "$RS" scripts/figures/progression_driver_coloc_phenotype_class_heatmap.R \
  > "$FIG2_SUPP_OUT_DIR/FigS2N.log" 2>&1; then
  exit 42
fi
grep -q 'FigS2N retired' "$FIG2_SUPP_OUT_DIR/FigS2N.log"
"$RS" scripts/figures/figS2O_pqtl_external_complementation.R \
  > "$FIG2_SUPP_OUT_DIR/FigS2O.log" 2>&1
"$RS" scripts/figures/figS2P_pqtl_eqtl_crossref.R \
  > "$FIG2_SUPP_OUT_DIR/FigS2P.log" 2>&1
"$RS" scripts/figures/figS2Q_scchromatin_external_corroboration.R \
  > "$FIG2_SUPP_OUT_DIR/FigS2Q.log" 2>&1
"$RS" scripts/figures/figS2R_mpra_external_corroboration.R \
  > "$FIG2_SUPP_OUT_DIR/FigS2R.log" 2>&1

find "$canonical_dir" -maxdepth 1 -type f -regextype posix-extended \
  -regex '.*/FigS2[J-R].*' -print0 | sort -z | xargs -0 -r sha256sum \
  > "$FIG2_SUPP_OUT_DIR/canonical_after.sha256"
cmp "$FIG2_SUPP_OUT_DIR/canonical_before.sha256" \
  "$FIG2_SUPP_OUT_DIR/canonical_after.sha256"

for f in "$FIG2_SUPP_OUT_DIR"/*.pdf; do
  pdfinfo "$f"
done > "$FIG2_SUPP_OUT_DIR/pdfinfo_all.txt"
if command -v pdffonts >/dev/null 2>&1; then
  for f in "$FIG2_SUPP_OUT_DIR"/*.pdf; do
    pdffonts "$f"
  done > "$FIG2_SUPP_OUT_DIR/pdffonts_all.txt"
else
  for f in "$FIG2_SUPP_OUT_DIR"/*.pdf; do
    grep -aEq '/FontFile[23]?' "$f"
    printf '%s\tembedded-font-object-present\n' "$(basename "$f")"
  done > "$FIG2_SUPP_OUT_DIR/pdffonts_all.txt"
fi
for f in "$FIG2_SUPP_OUT_DIR"/*.pdf; do
  pdftotext "$f" -
done > "$FIG2_SUPP_OUT_DIR/artwork_text.txt"

if grep -Eqi \
  'ABF-fallback|exploratory|causal gene|ancestry-specific|Cas13|expression-miss|maps null' \
  "$FIG2_SUPP_OUT_DIR/artwork_text.txt"; then
  exit 43
fi

pdf_count=$(find "$FIG2_SUPP_OUT_DIR" -maxdepth 1 -type f -name '*.pdf' | wc -l)
test "$pdf_count" -eq 7
page_count=$(grep -c '^Pages:[[:space:]]*1$' "$FIG2_SUPP_OUT_DIR/pdfinfo_all.txt")
test "$page_count" -eq 7
if command -v pdffonts >/dev/null 2>&1; then
  if awk 'NR > 2 && $0 ~ / no / { exit 1 }' "$FIG2_SUPP_OUT_DIR/pdffonts_all.txt"; then
    :
  else
    exit 44
  fi
fi

"$RS" scripts/figures/validate_figS2_JR_promoted_candidate.R \
  > "$FIG2_SUPP_OUT_DIR/validation.log" 2>&1
python3 scripts/manuscript/validate_resource_scope.py \
  > "$FIG2_SUPP_OUT_DIR/resource_scope.log" 2>&1
printf '%s\n' PASS > "$FIG2_SUPP_OUT_DIR/VALIDATION_STATUS.txt"
