#!/bin/bash
#SBATCH --job-name=format_panukbb
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/PanUKBB/format_panukbb_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/PanUKBB/format_panukbb_%j.err

# Format Pan-UKBB multi-ancestry GWAS (hg19→hg38 liftover)
# Also formats GLGC TG and MAGIC HOMA-IR for S-PrediXcan

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== Formatting Pipeline ==="
echo "Started: $(date)"
echo ""

# 1. Pan-UKBB: Extract AFR + CSA for ALT/AST/GGT and liftover
echo "--- Pan-UKBB ---"
cd "$BASE/GWAS/MR_Data/PanUKBB"
Rscript format_panukbb_for_coloc.R 2>&1

echo ""
echo "--- GLGC TG ---"
cd "$BASE/GWAS/MR_Data/GLGC"
Rscript format_glgc_for_spredixcan.R 2>&1

echo ""
echo "--- MAGIC HOMA-IR ---"
cd "$BASE/GWAS/MR_Data/MAGIC"
Rscript format_magic_for_spredixcan.R 2>&1

echo ""
echo "=== All formatting complete ==="
echo "Finished: $(date)"
