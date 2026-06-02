#!/bin/bash -l
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu
#SBATCH --array=1-22
#SBATCH --job-name=backfill_eur_bed
#SBATCH --output=GWAS/finemapping/logs/backfill_eur_bed_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/backfill_eur_bed_%A_%a.err

# Backfill .bed / .fam per 1kg_eur block. build_eur_ld.sh used --make-just-bim which
# produced only .ld + .bim. SuSiEX requires .bed/.fam too (shutil.copy2 each).
# We extract per-block PLINK BED using plink --make-bed on the chr-level file.

set -o pipefail
module load PLINK/1.9

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
cd "${BASE}"
CHR=${SLURM_ARRAY_TASK_ID}

EUR_DIR=data/ld_ref/1kg_eur
PLINK_PREFIX="${EUR_DIR}/chr${CHR}_eur"
BLOCKS_FILE="${EUR_DIR}/approx_LD_blocks.txt"

echo "[backfill] chr${CHR} start $(date)"
N_BLOCKS=0; N_SKIP=0; N_OK=0; N_FAIL=0

while read -r bchr bs be; do
  [ "${bchr}" = "chr" ] && continue
  [ "${bchr}" != "${CHR}" ] && continue
  N_BLOCKS=$((N_BLOCKS+1))
  BLOCK_DIR="${EUR_DIR}/chr${CHR}/${bs}.${be}"
  BLOCK_PREFIX="${BLOCK_DIR}/${bs}.${be}"

  # Skip if .bed already present
  if [ -f "${BLOCK_PREFIX}.bed" ] && [ -f "${BLOCK_PREFIX}.fam" ]; then
    N_SKIP=$((N_SKIP+1))
    continue
  fi

  # Needs: existing .bim (positions/variants to keep)
  if [ ! -f "${BLOCK_PREFIX}.bim" ]; then
    echo "  WARNING: no .bim for ${bs}.${be} — skipping"
    N_FAIL=$((N_FAIL+1))
    continue
  fi

  if plink --bfile "${PLINK_PREFIX}" \
       --chr "${CHR}" \
       --from-bp "${bs}" \
       --to-bp "${be}" \
       --make-bed \
       --out "${BLOCK_PREFIX}" \
       --allow-extra-chr \
       2>/dev/null; then
    [ -f "${BLOCK_PREFIX}.bed" ] && [ -f "${BLOCK_PREFIX}.fam" ] && N_OK=$((N_OK+1)) || N_FAIL=$((N_FAIL+1))
  else
    N_FAIL=$((N_FAIL+1))
  fi
done < "${BLOCKS_FILE}"

echo "[backfill] chr${CHR} done: ${N_OK} built, ${N_SKIP} pre-existing, ${N_FAIL} failed ($(date))"
