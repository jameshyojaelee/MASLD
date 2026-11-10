#!/bin/bash
# SS job S1d (per library, after STAR): IGH constant-region fragment counts for H38 (s3e_composition.R --ig).
# Reads only the IGH constant-exon span (chr14 and KI270846.1) of the sorted, indexed STAR BAM; the counting rule
# is the unchanged IGH block of extract_reads_one_library.py. Needs <bam>.bai and the STAR Log.final.out beside it.
# Writes <outdir>/<S>/{ig_fragments.tsv,library_qc.tsv,run_info.json}; refuses if <outdir>/<S> exists.
# Usage: s1d_ig_one.sh <S> <bam> <outdir>
set -euo pipefail
S=${1:?sample}; BAM=${2:?bam}; OUTDIR=${3:?outdir}
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
EXONS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/READS-20261009T141017Z/resources/final/ig_constant_exons.tsv
EXONS_SHA=08745694e2bd8d06de3c20c725cb2c4914df237692ce4a6f7d2b3f0ad6120e02
LOGF=${BAM%.Aligned.sortedByCoord.out.bam}.Log.final.out
[ -e "$OUTDIR/$S" ] && { echo "output exists; refusing: $OUTDIR/$S" >&2; exit 1; }
[ -s "$BAM.bai" ] || { echo "no index: $BAM.bai" >&2; exit 1; }
[ -s "$LOGF" ] || { echo "no STAR Log.final.out: $LOGF" >&2; exit 1; }
[ "$(sha256sum "$EXONS" | cut -d' ' -f1)" = "$EXONS_SHA" ] || { echo "ig_constant_exons.tsv sha256 mismatch" >&2; exit 1; }
mkdir -p "$OUTDIR"
"$PY" -I "$SS/scripts/s1d_ig_extract.py" --library "$S" --bam "$BAM" --log-final "$LOGF" --exons "$EXONS" --out-dir "$OUTDIR/$S"
