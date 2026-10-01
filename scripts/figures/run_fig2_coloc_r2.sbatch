#!/bin/bash
#SBATCH --job-name=fig2_coloc_r2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig2_coloc_r2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig2_coloc_r2_%j.err

# Figure 2 A-I against a corrected COLOC release (review item 2, Phase 3).
#
# Same steps as the canonical chain, which cannot be edited in place because it
# lives inside its own candidate directories:
#   figures/candidates/fig2_relayout_20260817T163114Z/render_fig2_relayout.sbatch
#     (variant classes + non-EUR GWS audit, rebuilt from the release master)
#   figures/candidates/fig2_coloc_audit_20260817T202451Z/render_and_validate_fig2.sbatch
#     (panels, page boxes, worked-locus and font checks, validator)
# The adopted-release literals (462-gene membership, RORA/FABP1 PP.H4, the sha256
# of the adopted master) are replaced by values read from the release's
# r2_summary.json and output_sha256.tsv.
#
#   sbatch --dependency=afterok:<coloc_r2 job> --export=ALL,\
# GENETICS_RELEASE_ROOT=<FM>/results/susie_coloc_r2_<ts>,\
# FIG2_R2_CANDIDATE=<BASE>/figures/candidates/fig2-coloc-r2-<UTC ts> \
#     scripts/figures/run_fig2_coloc_r2.sbatch
#
# Figure 2D reads pre-COLOC credible-set statistics; point FIG2_TRAIT_DIRECTNESS_STATS
# at the corrected review items 3-4 output when it exists (default: adopted).
set -euo pipefail
: "${GENETICS_RELEASE_ROOT:?Set the corrected COLOC release (results/susie_coloc_r2_<ts>)}"
: "${FIG2_R2_CANDIDATE:?Set a new figures/candidates/fig2-coloc-r2-<UTC ts> directory}"

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
CANDIDATE=$FIG2_R2_CANDIDATE
VARIANT_CLASS=$CANDIDATE/variant_class
COLOC=$GENETICS_RELEASE_ROOT/susie_coloc_all_gwas.csv
SUMMARY=$GENETICS_RELEASE_ROOT/r2_summary.json
R_RNASEQ=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
R_PLOT=/gpfs/commons/home/jameslee/micromamba/envs/plotgardener/bin/Rscript
PY_RNASEQ=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
RENDERER=/gpfs/commons/home/jameslee/.local/bin/pdftoppm
export PATH=/gpfs/commons/home/jameslee/.local/bin:$PATH

test ! -e "$CANDIDATE"
test -s "$GENETICS_RELEASE_ROOT/ASSEMBLY_COMPLETE.json"
test -s "$SUMMARY"
sha256sum -c "$GENETICS_RELEASE_ROOT/output_sha256.tsv"
mkdir -p "$VARIANT_CLASS"

cd "$BASE"
export MASLD_PROJECT_ROOT="$BASE"
export MASLD_FIGURE_PYTHON="$PY_RNASEQ"
export FIG2_CANDIDATE_DIR="$CANDIDATE"
export FIG2_VARIANT_CLASS_DIR="$VARIANT_CLASS"
export FIG2_COLOC_INPUT="$COLOC"
export FIG2_COLOC_SUMMARY="$SUMMARY"
export PYTHONNOUSERSITE=1

# Variant classes and the non-EUR GWS audit follow the release's positive set.
"$R_RNASEQ" scripts/figures/coloc_coding_noncoding/01_gather_variants.R
"$R_RNASEQ" scripts/figures/coloc_coding_noncoding/02_annotate_variants.R
"$R_RNASEQ" scripts/figures/coloc_coding_noncoding/03_aggregate_and_deg_overlap.R
"$R_PLOT" scripts/figures/audit_noneur_gws.R
export FIG2_NONEUR_AUDIT="$VARIANT_CLASS/noneur_gws_audit.csv"

FIG2_PANEL=cascade "$PY_RNASEQ" scripts/figures/fig1_gwas_creative_options.py \
  | tee "$CANDIDATE/fig2A_render.log"
grep '^CAPTION (Fig2A cascade)' "$CANDIDATE/fig2A_render.log" > "$CANDIDATE/fig2A_gwas_cascade_caption.txt"
"$R_RNASEQ" scripts/figures/fig2_resolution_plane.R
"$R_RNASEQ" scripts/figures/fig2_coding_noncoding.R

export FIG2_TRAIT_DIRECTNESS_STATS="${FIG2_TRAIT_DIRECTNESS_STATS:-$BASE/figures/candidates/fig2_trait_directness_panels_20260814T013923Z_codex}"
export FIG2_TRAIT_DIRECTNESS_OUT="$CANDIDATE"
"$R_PLOT" scripts/figures/fig2_trait_directness_panels.R

GATE=gws "$PY_RNASEQ" scripts/figures/fig2_ancestry_unique_coloc_gated.py
unset FIG2_CROSSANCESTRY_SOURCE
"$R_RNASEQ" scripts/figures/fig2_tri_ancestry_coloc.R

export LZ_MATCHED_MAIN=1
"$R_PLOT" scripts/figures/figS09_locus_zoom_disease.R \
  RORA "$CANDIDATE/Fig2G_RORA_NAFLD.pdf"
"$R_PLOT" scripts/figures/figS09_locus_zoom_disease.R \
  FABP1 "$CANDIDATE/Fig2H_FABP1_ALT.pdf"
"$R_PLOT" scripts/figures/fig2_locus_legend_strip.R \
  "$CANDIDATE/Fig2I_locus_legend.pdf"

"$R_RNASEQ" -e 'sessionInfo()' > "$CANDIDATE/sessionInfo_rnaseq.txt"
"$R_PLOT" -e 'sessionInfo()' > "$CANDIDATE/sessionInfo_plotgardener.txt"
"$PY_RNASEQ" -m pip freeze > "$CANDIDATE/pip_freeze.txt"
sha256sum "$COLOC" "$SUMMARY" "$VARIANT_CLASS/variant_classification.csv" \
  "$VARIANT_CLASS/noneur_gws_audit.csv" > "$CANDIDATE/input_sha256.tsv"
printf 'FIG2_TRAIT_DIRECTNESS_STATS\t%s\n' "$FIG2_TRAIT_DIRECTNESS_STATS" >> "$CANDIDATE/input_sha256.tsv"

"$PY_RNASEQ" - "$CANDIDATE" > "$CANDIDATE/page_box_audit.tsv" <<'PY'
import sys
from pathlib import Path
from pypdf import PdfReader

base = Path(sys.argv[1])
expected = {
    "fig2A_gwas_cascade.pdf": (3.53, 2.35),
    "Fig2B_PIP_vs_SuSiE-coloc.pdf": (1.80, 2.35),
    "Fig2C_coding_noncoding_split.pdf": (1.10, 2.35),
    "Fig2D_pip_architecture_by_trait_directness.pdf": (1.85, 2.10),
    "Fig2E_ancestry_unique_coloc_GWS.pdf": (2.00, 2.10),
    "Fig2F_crossancestry_coloc.pdf": (2.58, 2.10),
    "Fig2G_RORA_NAFLD.pdf": (3.26, 3.10),
    "Fig2H_FABP1_ALT.pdf": (3.26, 3.10),
    "Fig2I_locus_legend.pdf": (6.61, 0.28),
}
print("pdf\tpages\twidth_in\theight_in\texpected_w\texpected_h\tstatus")
for name, (want_w, want_h) in expected.items():
    reader = PdfReader(base / name)
    box = reader.pages[0].mediabox
    width, height = float(box.width) / 72, float(box.height) / 72
    passed = (
        len(reader.pages) == 1 and abs(width - want_w) < 1e-9
        and abs(height - want_h) < 1e-9
    )
    print(
        f"{name}\t{len(reader.pages)}\t{width:.3f}\t{height:.3f}\t"
        f"{want_w:.2f}\t{want_h:.2f}\t{'PASS' if passed else 'FAIL'}"
    )
    assert passed, name
PY

# Membership, untestable-gene caption and worked loci against the release.
"$PY_RNASEQ" - "$CANDIDATE" "$COLOC" "$SUMMARY" <<'PY'
import csv
import json
import math
import sys
from pathlib import Path

base, coloc = Path(sys.argv[1]), Path(sys.argv[2]).resolve()
summary = json.loads(Path(sys.argv[3]).read_text())

with (base / "Fig2E_ancestry_unique_coloc_GWS_membership.tsv").open() as handle:
    membership = list(csv.DictReader(handle, delimiter="\t"))
assert len(membership) == summary["fig2"]["multi_signal"], len(membership)

# Fig 2A: single-signal-only genes whose SuSiE test was untestable are listed in
# the source table and counted in the caption, with n from r2_summary.json.
with (base / "fig2A_gwas_cascade_source.tsv").open() as handle:
    cascade = [r for r in csv.DictReader(handle, delimiter="\t") if r["record_type"] == "union_gene"]
n_untestable = summary["fig2_susie_untestable"]["drawn_as_single_signal_only"]
assert len(cascade) == summary["fig2"]["multi_or_single_signal_coloc_gene_union"]
assert sum(r["signal_model"] == "single_signal_only" and r["susie_state"] == "untestable"
           for r in cascade) == n_untestable
caption = (base / "fig2A_gwas_cascade_caption.txt").read_text()
sentence = (f"{n_untestable} genes shown as single-signal only could not be tested by "
            "multi-signal colocalization")
assert (sentence in caption) == (n_untestable > 0), caption[-300:]

# The shared-variant PP is the ABF top_snp_PP, which the release does not change.
expected = {
    "Fig2G_RORA_NAFLD_source.csv": ("RORA", "MVP_NAFLD_EUR", 0.3754),
    "Fig2H_FABP1_ALT_source.csv": ("FABP1", "UKBB_ALT", 0.5578),
}
for name, (gene, gwas, shared_pp) in expected.items():
    locus = summary["worked_loci"][gene]
    assert locus["gwas_name"] == gwas and locus["method"] == "susie", locus
    with (base / name).open() as handle:
        row = next(csv.DictReader(handle))
    assert row["gene"] == gene and row["gwas"] == gwas
    assert math.isclose(float(row["coloc_PP4_susie"]), locus["PP.H4.susie"], abs_tol=5e-5)
    assert math.isclose(float(row["coloc_shared_variant_PP"]), shared_pp, abs_tol=5e-5)
    assert Path(row["coloc_input"]).resolve() == coloc
    assert row["coloc_method"] == "susie"
    assert row["gwas_ancestry"] == "EUR"
print("membership and worked loci match", sys.argv[3])
PY

for pdf in "$CANDIDATE"/*.pdf; do
  grep -aEq '/FontFile2?|/FontFile3' "$pdf"
  grep -aiEq 'Helvetica' "$pdf"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "$pdf"
done

# Printed values and rendered images for review of every panel.
mkdir -p "$CANDIDATE/previews"
for pdf in "$CANDIDATE"/*.pdf; do
  name=$(basename "$pdf" .pdf)
  pdftotext "$pdf" "$CANDIDATE/previews/$name.txt"
  "$RENDERER" -png -r 300 -singlefile "$pdf" "$CANDIDATE/previews/$name"
done

"$PY_RNASEQ" scripts/figures/validate_fig2_main_sequence.py \
  --panel-dir "$CANDIDATE" --layout-only --expected-json "$SUMMARY" \
  > "$CANDIDATE/figure2_validation.txt"
python3 scripts/manuscript/validate_resource_scope.py \
  > "$CANDIDATE/resource_scope_validation.txt"

echo "FIG2_COLOC_R2_COMPLETE $CANDIDATE"
