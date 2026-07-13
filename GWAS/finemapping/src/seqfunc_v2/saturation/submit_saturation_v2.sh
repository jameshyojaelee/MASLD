#!/usr/bin/env bash
# Submit only after the independent upstream gate can pass. The script never
# uses afterany: a failed gate or preparation prevents all expensive scoring.
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SRC=$ROOT/GWAS/finemapping/src/seqfunc_v2/saturation
OUT=${SATURATION_V2_ROOT:-$ROOT/GWAS/finemapping/results/seqfunc/haplotype_saturation/v2}
MODE=${1:-gate}
case "$MODE" in
  gate)
    sbatch "$SRC/00_gate.sbatch"
    ;;
  prepare)
    python3 "$SRC/00_upstream_gate.py"
    sbatch "$SRC/01_prepare.sbatch"
    ;;
  score)
    python3 "$SRC/00_upstream_gate.py"
    test -s "$OUT/sequences/scoring_shards.tsv" || { echo "run prepare first"; exit 2; }
    N=$(($(wc -l < "$OUT/sequences/scoring_shards.tsv")-1))
    test "$N" -gt 0
    SCORE=$(sbatch --parsable --array="0-$((N-1))%64" "$SRC/03_score.sbatch")
    sbatch --dependency="afterok:$SCORE" "$SRC/04_aggregate.sbatch"
    ;;
  *) echo "usage: $0 {gate|prepare|score}" >&2; exit 2;;
esac
