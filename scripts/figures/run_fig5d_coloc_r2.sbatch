#!/usr/bin/env bash
#SBATCH --job-name=fig5d_coloc_r2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig5d_coloc_r2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig5d_coloc_r2_%j.err

# Figure 5D on COLOC-eligible signal pairs only (review item 2, Phase 3).
# 1 restrict  sealed ATAC v3 context -> ATAC_R2_ROOT/genetics/context/ (primary
#             pair re-chosen among eligible pairs; no replay, no re-aggregation)
# 2 render    fig5d_variant_consistent_accessibility.R -> FIGURE_CANDIDATE_ROOT
# 3 check     panel counts rederived from the restricted table; printed n; PNG
#
#   sbatch --dependency=afterok:<coloc_r2 job> --export=ALL,\
# GENETICS_RELEASE_ROOT=<FM>/results/susie_coloc_r2_<ts>,\
# ATAC_R2_ROOT=<BASE>/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-<UTC ts>,\
# FIGURE_CANDIDATE_ROOT=<BASE>/figures/candidates/fig5d-coloc-r2-<UTC ts> \
#     scripts/figures/run_fig5d_coloc_r2.sbatch
#
# The full Figure 5 candidate (run_fig5_completion_candidate.sbatch and
# validate_fig5_completion.py) picks this up when FIG5D_PRIMARY_PAIRS and the
# three FIG5D_EXPECTED_* values printed at the end are exported.
set -euo pipefail
: "${GENETICS_RELEASE_ROOT:?Set the corrected COLOC release (results/susie_coloc_r2_<ts>)}"
: "${ATAC_R2_ROOT:?Set a new candidates/atac-context-v3-coloc-r2-<UTC ts> directory}"
: "${FIGURE_CANDIDATE_ROOT:?Set a new figures/candidates/fig5d-coloc-r2-<UTC ts> directory}"

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PATH=/gpfs/commons/home/jameslee/.local/bin:$PATH
test ! -e "$ATAC_R2_ROOT"
test ! -e "$FIGURE_CANDIDATE_ROOT"
test -s "$GENETICS_RELEASE_ROOT/ASSEMBLY_COMPLETE.json"
sha256sum -c "$GENETICS_RELEASE_ROOT/output_sha256.tsv"

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

"$PY_BIN" Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/24_restrict_primary_pairs_to_eligible.py \
  --out-root "$ATAC_R2_ROOT" --coloc-release "$GENETICS_RELEASE_ROOT"
CONTEXT="$ATAC_R2_ROOT/genetics/context"
census() { awk -F'\t' -v m="$1" '$1 == m {print $3}' "$CONTEXT/restriction_summary.tsv"; }
export FIG5D_PRIMARY_PAIRS="$CONTEXT/genetic_lineage_context_primary_pairs.tsv"
export FIG5D_EXPECTED_ROWS=$(census lineage_rows)
export FIG5D_EXPECTED_PAIRS=$(census gene_study_pairs)
export FIG5D_EXPECTED_GENES=$(census genes)
export FIGURE_CANDIDATE_ROOT

"$R_BIN" scripts/figures/fig5d_variant_consistent_accessibility.R

PANELS="$FIGURE_CANDIDATE_ROOT/figure5/panels"
PDF="$PANELS/fig5d_snatac_accessibility.pdf"
test "$(pdfinfo "$PDF" | awk '/^Pages:/ {print $2}')" = "1"
grep -aEq '/FontFile2?|/FontFile3' "$PDF"
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "$PDF"
mkdir -p "$FIGURE_CANDIDATE_ROOT/proofs"
pdftotext "$PDF" "$FIGURE_CANDIDATE_ROOT/proofs/fig5d_snatac_accessibility.txt"
pdftoppm -png -r 400 -singlefile "$PDF" "$FIGURE_CANDIDATE_ROOT/proofs/fig5d_snatac_accessibility"

# Every panel cell rederived from the restricted table; facet n printed on the panel.
"$PY_BIN" - "$FIG5D_PRIMARY_PAIRS" "$PANELS/data/fig5d_variant_consistent_accessibility.tsv" \
  "$FIGURE_CANDIDATE_ROOT/proofs/fig5d_snatac_accessibility.txt" <<'PY'
import csv, math, sys
from collections import Counter
primary_path, source_path, text_path = sys.argv[1:]
read = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
primary, source = read(primary_path), read(source_path)
observed = Counter((r["trait_class"], r["lineage"], r["evidence_state"]) for r in primary)
denominator = Counter(r["trait_class"] for r in primary if r["lineage"] == "hepatocyte")
assert len(source) == 50
for row in source:
    key = (row["trait_class"], row["lineage"], row["evidence_state"])
    assert int(row["n_pairs"]) == observed[key], key
    assert int(row["denominator"]) == denominator[row["trait_class"]], key
    assert math.isclose(float(row["proportion"]), int(row["n_pairs"]) / int(row["denominator"]), abs_tol=1e-12)
text = open(text_path).read()
for label, trait in (("MASLD/liver fat", "direct_MASLD"), ("Liver enzyme", "liver_enzyme")):
    assert f"{label} (n = {denominator[trait]:,})" in text, (label, denominator[trait])
print("fig5d cells rederived:", dict(denominator))
PY

sha256sum "$FIG5D_PRIMARY_PAIRS" "$PDF" "$PANELS/data/fig5d_variant_consistent_accessibility.tsv" \
  > "$FIGURE_CANDIDATE_ROOT/input_output_sha256.tsv"
FIG5_SESSION_INFO="$FIGURE_CANDIDATE_ROOT/sessionInfo.txt" \
  "$R_BIN" -e 'writeLines(capture.output(sessionInfo()), Sys.getenv("FIG5_SESSION_INFO"))'
python3 scripts/manuscript/validate_resource_scope.py > "$FIGURE_CANDIDATE_ROOT/resource_scope_validation.txt"
printf 'FIG5D_PRIMARY_PAIRS=%s FIG5D_EXPECTED_ROWS=%s FIG5D_EXPECTED_PAIRS=%s FIG5D_EXPECTED_GENES=%s\n' \
  "$FIG5D_PRIMARY_PAIRS" "$FIG5D_EXPECTED_ROWS" "$FIG5D_EXPECTED_PAIRS" "$FIG5D_EXPECTED_GENES"
echo "FIG5D_COLOC_R2_COMPLETE $FIGURE_CANDIDATE_ROOT"
