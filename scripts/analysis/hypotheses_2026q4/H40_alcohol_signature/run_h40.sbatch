#!/bin/bash
#SBATCH --job-name=H40alc
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
# H40: alcohol-associated signature in MASLD-labeled bulk biopsies.
# Usage: sbatch --output=<H40>/logs/slurm-%j.log run_h40.sbatch <H40_dir> <run_dir>
set -euo pipefail
H40="${1:?H40 dir required}"
OUT="${2:?run dir required}"
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H40_alcohol_signature
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER= PYTHONHASHSEED=0
test -s "$H40/PRESPEC_H40.md"
test ! -e "$OUT"
mkdir -p "$OUT"
cd "$OUT"
(cd "$SCRIPTS" && sha256sum parse_geo.py h40_alcohol_signature.R export_covariate_gene_sets.R run_h40.sbatch) > "$OUT/script_sha256.txt"
sha256sum "$H40"/downloads_untrusted/* "$H40"/gene_sets/*.tsv "$H40/PRESPEC_H40.md" > "$OUT/input_sha256.txt"
cp "$SCRIPTS/h40_alcohol_signature.R" "$OUT/h40_alcohol_signature.as_run.R"
cp "$SCRIPTS/parse_geo.py" "$OUT/parse_geo.as_run.py"
"$ENV/bin/python3" -I "$SCRIPTS/parse_geo.py" "$H40/downloads_untrusted" "$OUT"
"$ENV/bin/python3" -I -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
"$ENV/bin/python3" --version > "$OUT/python_version.txt"
"$ENV/bin/Rscript" --vanilla "$SCRIPTS/h40_alcohol_signature.R" "$H40" "$OUT" "${SLURM_CPUS_PER_TASK:-8}"
