#!/bin/bash
#SBATCH --job-name=fig3lnc
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3lnc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3lnc_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RNASEQ_ENV}/bin/Rscript" --vanilla scripts/figures/fig3_lncrna_role.R

PDF="${PROJECT_ROOT}/figures/main/fig3_bulk_transcriptomics/panels/fig3_lncrna_role_candidate.pdf"
test "$(pdfinfo "${PDF}" | awk '/^Pages:/ {print $2}')" = "1"
strings "${PDF}" | grep -q "/Subtype /TrueType"
! strings "${PDF}" | grep -q "/Subtype /Type3"
python3 scripts/manuscript/validate_resource_scope.py
