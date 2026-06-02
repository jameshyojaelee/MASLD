#!/usr/bin/env bash
# docs/build_html.sh — Render human-readable report docs from .md → .html
#
# Usage:
#   bash docs/build_html.sh           # rebuild all 6 report docs
#   bash docs/build_html.sh progress  # rebuild just progress.md
#
# Output: docs/html/<name>.html (standalone, embedded CSS, TOC, dark-mode aware)
# Source: docs/<name>.md (unchanged — .md stays the source of truth)

set -euo pipefail

PANDOC="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/.mamba/pathway_analysis/bin/pandoc"
CSS="docs/html/style.css"
OUTDIR="docs/html"

REPORT_DOCS=(
  paper_outline
  progress
  single_cell_analysis
  sex_stratified_analysis_framework
  dataset_labeling_and_harmonization
  SOP_Add_New_Dataset
)

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

if [[ $# -gt 0 ]]; then
  targets=("$@")
else
  targets=("${REPORT_DOCS[@]}")
fi

mkdir -p "$OUTDIR"

for name in "${targets[@]}"; do
  src="docs/${name}.md"
  out="${OUTDIR}/${name}.html"

  if [[ ! -f "$src" ]]; then
    echo "SKIP  $src (not found)" >&2
    continue
  fi

  "$PANDOC" "$src" \
    --standalone \
    --toc --toc-depth=3 \
    --number-sections \
    --css="style.css" \
    --mathjax \
    --metadata title="$(head -1 "$src" | sed 's/^#\+ *//')" \
    --metadata date="$(date -r "$src" '+%Y-%m-%d')" \
    --syntax-highlighting=tango \
    --wrap=none \
    -o "$out"

  size=$(du -h "$out" | cut -f1)
  echo "OK    $out  ($size)"
done
