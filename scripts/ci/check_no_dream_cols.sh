#!/usr/bin/env bash
# CI gate for the 2026-06-08 C2 canonical migration (dream -> limma-voom-qw).
# FAILS (exit 1) if any ACTIVE script still reads a retired dream_* atlas/consensus DEG
# column. The atlas now carries bulk_* (canonical_deg_results.csv uses unprefixed
# logFC/padj/t/lfsr). Intentional reads of the retired dream sensitivity arm must be
# annotated on the same line with the sentinel `# C2-OK-sensitivity` to pass.
#
# Benign / out-of-scope and therefore NOT flagged:
#   dream_robustness_flag (vestigial atlas flag), dream_pooled_results (mouse pipeline),
#   dream_comparator (proteomics column), dream_direction (Cas13 local var).
# Usage: bash scripts/ci/check_no_dream_cols.sh   (run from repo root)
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 2

PATTERN='dream_(padj|logFC_M|logFC_F|logFC|tstat|shrunk_logFC|lfsr|sig|dir)\b'

hits=$(grep -rEn \
  --include='*.R' --include='*.py' --include='*.sh' --include='*.js' --include='*.ts' --include='*.tsx' \
  --exclude-dir=archive --exclude-dir=.claude --exclude-dir=worktrees --exclude-dir=.git \
  --exclude-dir=.next --exclude-dir=out --exclude-dir=node_modules --exclude-dir=dist \
  --exclude='check_no_dream_cols.sh' \
  "$PATTERN" . 2>/dev/null \
  | grep -vE 'docs/audit/|docs/archive/|_backup|\.bak:|/\.next/|/out/|/node_modules/|RNA-seq/Mouse/' \
  | grep -v 'C2-OK-sensitivity' \
  || true)

if [ -n "$hits" ]; then
  n=$(printf '%s\n' "$hits" | wc -l | tr -d ' ')
  echo "FAIL: $n unannotated retired dream_* column reference(s) remain in active code:"
  printf '%s\n' "$hits"
  exit 1
fi
echo "PASS: no unannotated retired dream_* column references in active code."
exit 0
