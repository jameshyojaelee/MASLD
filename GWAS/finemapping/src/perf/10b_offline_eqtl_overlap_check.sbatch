#!/bin/bash
#SBATCH --job-name=coloc_elig_offline
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_eligibility/offline_check_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_eligibility/offline_check_%j.err

# Offline eQTL shared-posterior check of the 816 tier-1/2 SuSiE-positive pairs
# replayed by the ATAC v3 candidate (review item 2). No SuSiE refit.
#   sbatch 10b_offline_eqtl_overlap_check.sbatch <NEW_OUT_DIR>
set -eo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch 10b_offline_eqtl_overlap_check.sbatch <NEW_OUT_DIR>" >&2
  exit 2
fi
OUT_DIR=$1
if [[ -e "${OUT_DIR}" ]]; then echo "refusing to reuse ${OUT_DIR}" >&2; exit 1; fi
mkdir -p "${OUT_DIR}"

set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
export EQTL_SUSIE_DIR="${FM}/results/eqtl_susie_polyfun"
OUT_TSV="${OUT_DIR}/eqtl_overlap_check.tsv" Rscript "${FM}/src/perf/10_offline_eqtl_overlap_check.R"
