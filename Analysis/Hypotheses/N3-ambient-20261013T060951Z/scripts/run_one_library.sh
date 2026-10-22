#!/bin/bash
# One library: export to node-local temporary files, decontX, keep only the small result.
# Usage: run_one_library.sh <LIBRARY> <DATASET>
set -euo pipefail
LIB="$1"; DS="$2"
OUT="$RES/per_library/$LIB.rds"
if [ -f "$OUT" ]; then echo "[skip] $LIB exists"; exit 0; fi
"$PY" -I "$SCR/01_export_library.py" "$REPO" "$LIB" "$DS" "$WORK"
"$RS_DECONTX" "$SCR/02_decontx_library.R" "$WORK/$LIB" "$OUT" "$RES/report_genes.txt" "$RES/genes.txt"
rm -r "${WORK:?}/$LIB"
