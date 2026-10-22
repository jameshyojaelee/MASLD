#!/bin/bash
#SBATCH --job-name=H49rma
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
# H49 step B: per-compound RMA + limma on the TG-GATEs human in vitro archives.
# Usage: sbatch --output=<H49>/logs/rma-%j.log run_h49_rma.sbatch <H49_dir> <out_dir> [compound ...]
set -euo pipefail
H49="${1:?H49 dir required}"
OUT="${2:?out dir required}"
shift 2
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H49_drug_steatosis_early_genes
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER=
test -s "$H49/PRESPEC_H49.md"
test ! -e "$OUT"
mkdir -p "$OUT"
JOBTMP=$(mktemp -d "/tmp/h49_${SLURM_JOB_ID:-local}_XXXX")
export TMPDIR="$JOBTMP"
cd "$SCRIPTS"
sha256sum h49_rma_limma.R > "$OUT/script_sha256.txt"
"$ENV/bin/Rscript" --vanilla h49_rma_limma.R \
  "$H49/downloads_untrusted_tggates_cel" \
  "$H49/downloads_untrusted_tggates_meta_unzipped/open_tggates_cel_file_attribute.csv" \
  "$H49/downloads_untrusted_annotation/hgu133plus2cdf_2.18.0.tar.gz" \
  "$OUT" "$@"
rmdir "$JOBTMP" 2>/dev/null || true
