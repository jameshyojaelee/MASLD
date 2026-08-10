#!/bin/bash
#SBATCH --job-name=susie_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_gram/coloc_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_gram/coloc_%A_%a.err

# Re-run the canonical SuSiE-COLOC pipeline for the 29 NON-EUR GWAS against the
# REBUILT (mean-imputed Gram) 1000G LD panels.
#
# Why: get_ld_per_locus() reads the same block .ld files the seqfunc work showed
# to be indefinite, and 06_susie_coloc.R adds a fixed LD_REGULARIZE = 1e-3.
# Measured on an AFR block, chol(R + 1e-3 I) succeeds for 78% of 200-variant
# cis-windows but 0% at >=1000 variants -- so the SuSiE arm silently fails into
# its tryCatch and the gene falls back to ABF-only. The rebuilt panels succeed
# 100% at every window size tested.
#
# SAFETY: COLOC_OUT_SUFFIX=_gram writes to results/susie_coloc_gram/<gwas>/,
# side by side with the canonical results/susie_coloc/. NOTHING canonical is
# overwritten -- gene_level_coloc.csv is firewall-guarded and feeds the atlas,
# so promotion is a separate, explicit decision.
#
# The ABF arm is unaffected by any of this (it uses no LD) and is recomputed
# here only because it shares the same script; it should reproduce.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SRC="${PROJECT_ROOT}/GWAS/finemapping/src"
LDREF="${PROJECT_ROOT}/GWAS/finemapping/data/ld_ref"
export MASLD_PROJECT_ROOT

# Per-ancestry overrides beat the LD_PANEL dispatch in get_ld_base_dir().
export AFR_LD_DIR="${LDREF}/1kg_afr_gram"
export AMR_LD_DIR="${LDREF}/1kg_amr_gram"
export EAS_LD_DIR="${LDREF}/1kg_eas_gram"
export SAS_LD_DIR="${LDREF}/1kg_sas_gram"
export COLOC_OUT_SUFFIX="_gram"

export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

GWAS_LIST=(
  BBJ_ALT BBJ_AST BBJ_GGT
  MVP_Albumin_AFR MVP_Albumin_AMR MVP_Albumin_EAS
  MVP_ALT_AFR MVP_ALT_AMR MVP_ALT_EAS
  MVP_AST_AFR MVP_AST_AMR MVP_AST_EAS
  MVP_ChronLiver_AFR MVP_ChronLiver_AMR MVP_ChronLiver_EAS
  MVP_Cirrhosis_AFR MVP_Cirrhosis_AMR
  MVP_NAFLD_AFR MVP_NAFLD_AMR MVP_NAFLD_EAS
  MVP_Platelet_AFR MVP_Platelet_AMR MVP_Platelet_EAS
  PanUKBB_AFR_ALT PanUKBB_AFR_AST PanUKBB_AFR_GGT
  PanUKBB_CSA_ALT PanUKBB_CSA_AST PanUKBB_CSA_GGT
)
N_GWAS=${#GWAS_LIST[@]}

# Array index -> (gwas, chromosome).  Task IDs are 0-based over 29 x 22.
IDX=${SLURM_ARRAY_TASK_ID}
GI=$(( IDX / 22 ))
CHR=$(( IDX % 22 + 1 ))
GWAS=${GWAS_LIST[$GI]}

echo "[$(date -u +%FT%TZ)] task ${IDX} -> ${GWAS} chr${CHR}  (gram LD, out suffix _gram)"
cd "${SRC}"
GWAS_NAME="${GWAS}" CHR_FILTER="${CHR}" Rscript "${SRC}/06_susie_coloc.R"
echo "[$(date -u +%FT%TZ)] done ${GWAS} chr${CHR}"
