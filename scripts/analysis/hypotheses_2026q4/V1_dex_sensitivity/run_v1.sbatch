#!/bin/bash
#SBATCH --job-name=V1dex
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
# V1: GSE320593 dexamethasone response vs bulk disease signature.
# Usage: sbatch --output=<V1>/logs/slurm-%j.log run_v1.sbatch <V1_dir> <out_dir>
set -euo pipefail
V1="${1:?V1 dir required}"
OUT="${2:?out dir required}"
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/V1_dex_sensitivity
export OMP_NUM_THREADS=2 R_LIBS_USER= PYTHONHASHSEED=0
test -s "$V1/PRESPEC_V1.md"
test ! -e "$OUT"
mkdir -p "$OUT"
cd "$OUT"
(cd "$SCRIPTS" && sha256sum parse_soft.py v1_dex_overlap.R run_v1.sbatch) > "$OUT/script_sha256.txt"
sha256sum "$V1"/downloads_untrusted/GSE320593_* "$V1/PRESPEC_V1.md" > "$OUT/input_sha256.txt"
cp "$SCRIPTS/v1_dex_overlap.R" "$OUT/v1_dex_overlap.as_run.R"
cp "$SCRIPTS/parse_soft.py" "$OUT/parse_soft.as_run.py"
"$ENV/bin/python3" -I "$SCRIPTS/parse_soft.py" \
  "$V1/downloads_untrusted/GSE320593_family.soft.gz" "$OUT/sample_design.tsv"
"$ENV/bin/python3" --version > "$OUT/python_version.txt"
"$ENV/bin/Rscript" --vanilla "$SCRIPTS/v1_dex_overlap.R" "$V1" "$OUT" "$OUT/sample_design.tsv"
