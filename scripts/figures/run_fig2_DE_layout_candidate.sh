#!/usr/bin/env bash
set -euo pipefail

: "${MASLD_PROJECT_ROOT:?MASLD_PROJECT_ROOT is required}"
: "${FIG2_DE_OUT:?FIG2_DE_OUT is required}"
: "${FIG2_TRAIT_DIRECTNESS_STATS:?FIG2_TRAIT_DIRECTNESS_STATS is required}"

mkdir -p "$FIG2_DE_OUT"
cd "$MASLD_PROJECT_ROOT"
printf 'job_id\tstarted_utc\n%s\t%s\n' \
  "${SLURM_JOB_ID:-not_slurm}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  > "$FIG2_DE_OUT/PREFLIGHT_STARTED.tsv"

canonical_dir="figures/main/fig2_genetics/panels"
sha256sum "$canonical_dir/Fig2D_pip_architecture_by_trait_directness.pdf" \
  "$canonical_dir/Fig2E_ancestry_unique_coloc_GWS.pdf" \
  > "$FIG2_DE_OUT/canonical_before.sha256"

export FIG2_TRAIT_DIRECTNESS_OUT="$FIG2_DE_OUT"
export FIG2_CANDIDATE_DIR="$FIG2_DE_OUT"
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq_r44/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export MASLD_FIGURE_PYTHON="$PY"
export PYTHONNOUSERSITE=1
export FIG2_COLOC_INPUT="$MASLD_PROJECT_ROOT/GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"
export FIG2_NONEUR_AUDIT="$MASLD_PROJECT_ROOT/figures/candidates/fig2_relayout_20260817T164455Z/variant_class/noneur_gws_audit.csv"
test "$(sha256sum "$FIG2_COLOC_INPUT" | cut -d' ' -f1)" = \
  a78ec6bd7ffd93a49b5c6bfa291eb3fb9ca3dc84a9e88b24878a2bf44ecfcbc9
test "$(sha256sum "$FIG2_NONEUR_AUDIT" | cut -d' ' -f1)" = \
  d132762c5b7868626524208df7d227f5fe5a00d4ac9f67aa91f03201306dff5c
"$R" scripts/figures/fig2_trait_directness_panels.R \
  > "$FIG2_DE_OUT/Fig2D_render.log" 2>&1
GATE=gws "$PY" scripts/figures/fig2_ancestry_unique_coloc_gated.py \
  > "$FIG2_DE_OUT/Fig2E_render.log" 2>&1

sha256sum "$canonical_dir/Fig2D_pip_architecture_by_trait_directness.pdf" \
  "$canonical_dir/Fig2E_ancestry_unique_coloc_GWS.pdf" \
  > "$FIG2_DE_OUT/canonical_after.sha256"
cmp "$FIG2_DE_OUT/canonical_before.sha256" "$FIG2_DE_OUT/canonical_after.sha256"

for pdf in \
  "$FIG2_DE_OUT/Fig2D_pip_architecture_by_trait_directness.pdf" \
  "$FIG2_DE_OUT/Fig2E_ancestry_unique_coloc_GWS.pdf"; do
  pdfinfo "$pdf"
  pdftotext "$pdf" -
  grep -aEq '/FontFile[23]?' "$pdf"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "$pdf"
done > "$FIG2_DE_OUT/pdf_validation.txt"

"$PY" - "$FIG2_DE_OUT" <<'PY'
import re
import subprocess
import sys
from pathlib import Path

out = Path(sys.argv[1])
expected = {
    "Fig2D_pip_architecture_by_trait_directness.pdf": (1.85, 2.10),
    "Fig2E_ancestry_unique_coloc_GWS.pdf": (2.00, 2.10),
}
for name, (want_w, want_h) in expected.items():
    pdf = out / name
    info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
    assert re.search(r"^Pages:\s+1$", info, re.M), info
    size = re.search(r"^Page size:\s+([\d.]+) x ([\d.]+) pts", info, re.M)
    assert size
    got_w, got_h = float(size.group(1)) / 72, float(size.group(2)) / 72
    assert abs(got_w - want_w) <= 0.002, (name, got_w, want_w)
    assert abs(got_h - want_h) <= 0.002, (name, got_h, want_h)

d = (out / "Fig2D_layout_measurements.tsv").read_text()
e = (out / "Fig2E_layout_measurements.tsv").read_text()
assert "right_artwork_clearance_in" in d
assert "x_label_horizontal_gap_in" in d
assert "footnote_to_small_label_vertical_gap_in" in e
print("PASS")
PY

python3 scripts/manuscript/validate_resource_scope.py \
  > "$FIG2_DE_OUT/resource_scope.log" 2>&1
printf '%s\n' PASS > "$FIG2_DE_OUT/VALIDATION_STATUS.txt"
