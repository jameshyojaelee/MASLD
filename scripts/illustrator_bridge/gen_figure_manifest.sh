#!/usr/bin/env bash
# gen_figure_manifest.sh
# RUN ON THE HPC. Writes figures/figures_manifest.tsv — a lightweight index
# (relative path, size, mtime) of every figure PDF, so the Mac side can:
#   (1) sync efficiently, and
#   (2) know which panels exist and when each last changed.
# Light enough for the login node (a single find).
set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
FIG_DIR="$BASE/figures"
OUT="$FIG_DIR/figures_manifest.tsv"

{
  printf 'rel_path\tbytes\tmtime_epoch\tmtime_human\n'
  find "$FIG_DIR" -name '*.pdf' -type f \
    -printf '%P\t%s\t%T@\t%TY-%Tm-%Td %TH:%TM\n' | sort
} > "$OUT"

echo "Wrote $OUT ($(( $(wc -l < "$OUT") - 1 )) PDFs)"
