#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --partition=cpu
#SBATCH --array=1-22
#SBATCH --job-name=build_eur_ld
#SBATCH --output=logs/build_eur_ld_%A_%a.out
#SBATCH --error=logs/build_eur_ld_%A_%a.err

# build_eur_ld.sh
# Builds EUR per-block LD matrices from 1000G Phase 1 PLINK files in the same
# block structure as the existing 1kg_{eas,afr,sas} layout (Berisa-Pickrell EUR
# blocks, 1704 total). One array task per chromosome.
#
# Inputs:
#   - GWAS/finemapping/data/ld_ref/1kg_eur/chr${CHR}_eur.{bed,bim,fam}
#       (symlinks to data/1kg_eur/bed/chr${CHR}_eur.* — built by RNA-seq/prep_1kg_eur_ld.sh)
#   - GWAS/finemapping/data/ld_ref/1kg_eur/approx_LD_blocks.txt
#       (Berisa-Pickrell EUR; installed 2026-04-23)
#
# Outputs:
#   - GWAS/finemapping/data/ld_ref/1kg_eur/chr${CHR}/{start}.{stop}/{start}.{stop}.ld
#   - GWAS/finemapping/data/ld_ref/1kg_eur/chr${CHR}/{start}.{stop}/{start}.{stop}.bim
#
# Note: 1kg_eur has 379 EUR samples (much smaller than UKBB's 337K). LD estimates
# will be noisier — this is the cost of severing sghatan dependency. TOP-LD EUR
# (n=13,160) provides a higher-quality alternative; see Phase 2.

set -euo pipefail

module load PLINK/1.9

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

CHR=${SLURM_ARRAY_TASK_ID}
echo "============================================================"
echo "Building EUR LD matrices for chr${CHR}"
echo "============================================================"

EUR_DIR="data/ld_ref/1kg_eur"
PLINK_PREFIX="${EUR_DIR}/chr${CHR}_eur"
BLOCKS_FILE="${EUR_DIR}/approx_LD_blocks.txt"

# Sanity checks
for ext in bed bim fam; do
  if [ ! -f "${PLINK_PREFIX}.${ext}" ]; then
    echo "ERROR: input PLINK file missing: ${PLINK_PREFIX}.${ext}"
    exit 1
  fi
done
if [ ! -f "${BLOCKS_FILE}" ]; then
  echo "ERROR: block file missing: ${BLOCKS_FILE}"
  exit 1
fi

N_VARS_TOTAL=$(wc -l < "${PLINK_PREFIX}.bim")
N_SAMPLES=$(wc -l < "${PLINK_PREFIX}.fam")
echo "  Input: ${N_VARS_TOTAL} variants, ${N_SAMPLES} EUR samples"

mkdir -p "${EUR_DIR}/chr${CHR}"

N_BLOCKS=0
N_SUCCESS=0
N_SKIPPED=0
N_FAILED=0

while read -r block_chr block_start block_stop; do
  # Skip header line and other chromosomes
  if [ "${block_chr}" != "${CHR}" ] || [ "${block_chr}" == "chr" ]; then
    continue
  fi

  N_BLOCKS=$((N_BLOCKS + 1))
  BLOCK_DIR="${EUR_DIR}/chr${CHR}/${block_start}.${block_stop}"
  BLOCK_PREFIX="${BLOCK_DIR}/${block_start}.${block_stop}"

  # Skip if already computed
  if [ -f "${BLOCK_PREFIX}.ld" ] && [ -f "${BLOCK_PREFIX}.bim" ]; then
    N_SKIPPED=$((N_SKIPPED + 1))
    N_SUCCESS=$((N_SUCCESS + 1))
    continue
  fi

  mkdir -p "${BLOCK_DIR}"

  # Compute pairwise R for this block
  if plink --bfile "${PLINK_PREFIX}" \
      --chr "${CHR}" \
      --from-bp "${block_start}" \
      --to-bp "${block_stop}" \
      --r square \
      --make-just-bim \
      --out "${BLOCK_PREFIX}" \
      --allow-extra-chr \
      2>/dev/null; then
    if [ -f "${BLOCK_PREFIX}.ld" ] && [ -f "${BLOCK_PREFIX}.bim" ]; then
      N_SUCCESS=$((N_SUCCESS + 1))
      N_VARS=$(wc -l < "${BLOCK_PREFIX}.bim")
      if [ $((N_BLOCKS % 10)) -eq 0 ]; then
        echo "  Block ${N_BLOCKS}: ${block_start}-${block_stop} (${N_VARS} variants)"
      fi
    else
      N_FAILED=$((N_FAILED + 1))
      echo "  WARNING: PLINK succeeded but output files missing for block ${block_start}-${block_stop}"
    fi
  else
    N_FAILED=$((N_FAILED + 1))
    echo "  WARNING: PLINK failed for block ${block_start}-${block_stop} (likely no variants in range)"
  fi
done < "${BLOCKS_FILE}"

echo ""
echo "chr${CHR}: ${N_SUCCESS}/${N_BLOCKS} blocks computed (${N_SKIPPED} skipped pre-existing, ${N_FAILED} failed)"
echo "============================================================"
echo "Done building 1kg_eur LD for chr${CHR}: $(date)"
echo "============================================================"
