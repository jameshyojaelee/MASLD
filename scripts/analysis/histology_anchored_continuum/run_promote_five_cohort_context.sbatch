#!/usr/bin/env bash
#SBATCH --job-name=hac_promote_5cohort
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PROMOTER="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/promote_five_cohort_context.py"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

python3 -m py_compile "${PROMOTER}"
python3 "${PROMOTER}"

for manifest in \
  "${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/manifests/continuum_five_cohort_revision_20260824.tsv" \
  "${PROJECT_ROOT}/figures/supplementary/figS03_bulk_transcriptomics/manifests/continuum_five_cohort_revision_20260824.tsv"; do
  while IFS=$'\t' read -r callout destination source old_sha new_sha rest; do
    [[ "${callout}" == "callout" ]] && continue
    pdf="${PROJECT_ROOT}/${destination}"
    [[ "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" == "1" ]]
    gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
    ! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${pdf}"
    [[ "$(sha256sum "${pdf}" | awk '{print $1}')" == "${new_sha}" ]]
  done < "${manifest}"
done

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
