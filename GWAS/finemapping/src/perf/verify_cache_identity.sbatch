#!/bin/bash
#SBATCH --job-name=susieR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/verify_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/verify_%A_%a.err
#
# IDENTITY + SPEED verification for the LD block cache (perf/get_ld_per_locus_cached.R).
#
# Target: 2021_34128465_PDFF_EUR chr21 -- array task 64 of the live rerun, which
# COMPLETED in 2.89h. Its output already exists, so both arms below have a fixed
# reference to be compared against.
#
#   task 0 = CONTROL : canonical 06_susie_coloc.R, unpatched  -> suffix _verifyctl
#   task 1 = CACHED  : same script + the override             -> suffix _verifycache
#
# WHY A CONTROL ARM: without it, any difference in the cached arm is ambiguous --
# it could be my bug OR ordinary nondeterminism in the pipeline. The control
# establishes whether re-running the SAME code reproduces the SAME bytes. Only if
# CONTROL == reference is a CACHED == reference comparison meaningful.
#
# Wall-clock of the two arms is the speedup measurement: identical work, one
# variable changed.
#
# TOUCHES NOTHING LIVE: writes to results/susie_coloc_verifyctl|_verifycache.
# The canonical finemapping_functions.R is not modified -- the override shadows
# it by being sourced afterwards, and only in the patched copy.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
FM="${PROJECT_ROOT}/GWAS/finemapping"
SRC="${FM}/src"
LDREF="${FM}/data/ld_ref"
export MASLD_PROJECT_ROOT

# Identical environment to 06c_coloc_rerun_full.sbatch, so the only variable is the code.
export AFR_LD_DIR="${LDREF}/1kg_afr_gram"
export AMR_LD_DIR="${LDREF}/1kg_amr_gram"
export EAS_LD_DIR="${LDREF}/1kg_eas_gram"
export SAS_LD_DIR="${LDREF}/1kg_sas_gram"
export EQTL_SUSIE_DIR="${FM}/results/eqtl_susie_polyfun"
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

GWAS="2021_34128465_PDFF_EUR"
CHR=21

if [ "${SLURM_ARRAY_TASK_ID}" = "0" ]; then
  ARM="control"; SCRIPT="${SRC}/06_susie_coloc.R";              export COLOC_OUT_SUFFIX="_verifyctl"
else
  ARM="cached";  SCRIPT="${SRC}/perf/06_susie_coloc_cached.R";  export COLOC_OUT_SUFFIX="_verifycache"
fi

echo "[$(date -u +%FT%TZ)] arm=${ARM} script=${SCRIPT} suffix=${COLOC_OUT_SUFFIX}"
echo "reference = ${FM}/results/susie_coloc_rerun/${GWAS}/susie_coloc_chr${CHR}.csv"

# HARD GUARD. The rm -f lines below interpolate COLOC_OUT_SUFFIX into a path. If
# that variable were ever empty the target would resolve to results/susie_coloc/
# -- the CANONICAL directory that feeds gene_level_coloc.csv, the atlas and Fig2.
# Refuse to run at all rather than rely on the assignment above staying correct.
case "${COLOC_OUT_SUFFIX}" in
  _verify*) : ;;
  *) echo "FATAL: COLOC_OUT_SUFFIX='${COLOC_OUT_SUFFIX}' is not a _verify* suffix; refusing to touch anything" >&2; exit 1 ;;
esac

# Fresh start: a stale checkpoint would let the arm skip the very genes under test.
rm -f "${FM}/results/susie_coloc${COLOC_OUT_SUFFIX}/${GWAS}/checkpoint_chr${CHR}.csv" 2>/dev/null || true
rm -f "${FM}/results/susie_coloc${COLOC_OUT_SUFFIX}/${GWAS}/susie_coloc_chr${CHR}.csv" 2>/dev/null || true

cd "${SRC}"
START=$(date +%s)
GWAS_NAME="${GWAS}" CHR_FILTER="${CHR}" Rscript "${SCRIPT}"
END=$(date +%s)

echo "[$(date -u +%FT%TZ)] arm=${ARM} done in $(( END - START )) s"

REF="${FM}/results/susie_coloc_rerun/${GWAS}/susie_coloc_chr${CHR}.csv"
OUT="${FM}/results/susie_coloc${COLOC_OUT_SUFFIX}/${GWAS}/susie_coloc_chr${CHR}.csv"
echo "=== DIFF arm=${ARM} vs reference ==="
if [ ! -f "$OUT" ]; then
  echo "RESULT ${ARM}: NO OUTPUT PRODUCED"
elif diff -q "$REF" "$OUT" >/dev/null 2>&1; then
  echo "RESULT ${ARM}: BYTE-IDENTICAL to reference"
else
  echo "RESULT ${ARM}: DIFFERS from reference"
  echo "  ref rows=$(($(wc -l < "$REF")-1))  out rows=$(($(wc -l < "$OUT")-1))"
  diff "$REF" "$OUT" | head -20
fi
echo "SECONDS ${ARM} $(( END - START ))"
