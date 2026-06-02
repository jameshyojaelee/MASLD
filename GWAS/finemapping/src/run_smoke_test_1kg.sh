#!/bin/bash -l
#SBATCH --job-name=smoke_1kg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=GWAS/finemapping/logs/smoke_1kg_%j.out
#SBATCH --error=GWAS/finemapping/logs/smoke_1kg_%j.err

# Smoke test: re-run UKBB_GGT chr15 SuSiE-COLOC with LD_PANEL=1kg.
# Compare RORA PP.H4.susie to v1 UKBB baseline (0.9954).

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${BASE}"

# Output via env-var redirect — writes to results/susie_coloc_1kg/UKBB_GGT/
export GWAS_NAME=UKBB_GGT
export CHR_FILTER=15
export LD_PANEL=1kg
export COLOC_OUT_SUFFIX="_1kg"

SMOKE_DIR="${BASE}/GWAS/finemapping/results/susie_coloc_1kg/UKBB_GGT"

echo "[smoke] start $(date)"
echo "[smoke] LD_PANEL=${LD_PANEL}  COLOC_OUT_SUFFIX=${COLOC_OUT_SUFFIX}"
echo "[smoke] running 06_susie_coloc.R UKBB_GGT 15"

Rscript GWAS/finemapping/src/06_susie_coloc.R UKBB_GGT 15

echo "[smoke] done $(date)"
echo ""
echo "=== RORA result (1kg LD) ==="
grep "^RORA" "${SMOKE_DIR}/susie_coloc_chr15.csv" 2>/dev/null
echo ""
echo "=== RORA result (v1 UKBB baseline for comparison) ==="
grep "^RORA" "${BASE}/GWAS/finemapping/results/susie_coloc_ukbb_sghatan_v1_2026-04-21/UKBB_GGT/susie_coloc_chr15.csv"
