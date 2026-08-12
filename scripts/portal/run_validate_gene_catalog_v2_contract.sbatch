#!/bin/bash
#SBATCH --job-name=catreview
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/catreview_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/catreview_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PYTHON="/gpfs/commons/home/jameslee/micromamba/envs/alphagenome/bin/python"
CANDIDATE_DIR="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/gene-catalog-v2-contract-2026-08-11"
REVIEW_DIR="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/gene-catalog-v2-contract-review-2026-08-11"

cd "${PROJECT_ROOT}"
PYTHONPATH="${PROJECT_ROOT}/scripts/portal" "${PYTHON}" \
  scripts/portal/validate_gene_catalog_v2_contract.py \
  --candidate-dir "${CANDIDATE_DIR}" \
  --review-dir "${REVIEW_DIR}"
