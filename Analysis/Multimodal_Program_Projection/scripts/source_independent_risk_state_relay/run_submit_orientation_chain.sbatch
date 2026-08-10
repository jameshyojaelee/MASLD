#!/bin/bash
#SBATCH --job-name=orientsub
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/orientsub_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/orientsub_%j.err

set -euo pipefail

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
CANDIDATE_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates
: "${PLAN45_ORIENTATION_CANDIDATE_ID:?PLAN45_ORIENTATION_CANDIDATE_ID is required}"
: "${PLAN45_ROUTING_CANDIDATE_ID:?PLAN45_ROUTING_CANDIDATE_ID is required}"
: "${PLAN45_GUIDE_CANDIDATE_ID:?PLAN45_GUIDE_CANDIDATE_ID is required}"

orientation_root=${CANDIDATE_ROOT}/${PLAN45_ORIENTATION_CANDIDATE_ID}
routing_root=${CANDIDATE_ROOT}/${PLAN45_ROUTING_CANDIDATE_ID}
guide_root=${CANDIDATE_ROOT}/${PLAN45_GUIDE_CANDIDATE_ID}
worklist=${orientation_root}/orientation_worklist.tsv
seal=${orientation_root}/ORIENTATION_WORKLIST_SEALED.json
lineage_root=${CANDIDATE_ROOT}/source-independent-risk-state-relay-lineage-observability-2026-08-10
accessibility_root=${CANDIDATE_ROOT}/source-independent-risk-state-relay-accessibility-reference-2026-08-10

if [[ ! -s "${worklist}" || ! -s "${seal}" ]]; then
  echo "sealed orientation worklist is absent"
  exit 1
fi
if [[ -e "${orientation_root}/ORIENTATION_RELEASE.json" ]]; then
  echo "orientation release already exists; refusing to resubmit"
  exit 1
fi
if [[ -e "${routing_root}" || -e "${guide_root}" ]]; then
  echo "downstream candidate identity already exists; refusing to overwrite"
  exit 1
fi
for dependency in \
  "${lineage_root}/LINEAGE_REFERENCE_SEALED.json" \
  "${accessibility_root}/ACCESSIBILITY_REFERENCE_SEALED.json"; do
  if [[ ! -s "${dependency}" ]]; then
    echo "required frozen routing dependency is absent: ${dependency}"
    exit 1
  fi
done

n_tasks=$(awk 'END {print NR - 1}' "${worklist}")
if [[ ! "${n_tasks}" =~ ^[1-9][0-9]*$ ]]; then
  echo "orientation worklist has no executable rows: ${n_tasks}"
  exit 1
fi
last_index=$((n_tasks - 1))
cd "${SCRIPT_ROOT}"

pair_job_raw=$(sbatch --parsable \
  --array="0-${last_index}%20" \
  --export="ALL,MASLD_PROJECT_ROOT=${PROJECT_ROOT},PLAN45_CANDIDATE_ID=${PLAN45_ORIENTATION_CANDIDATE_ID}" \
  "${SCRIPT_ROOT}/run_pair_export.sbatch")
pair_job=${pair_job_raw%%;*}
validator_job_raw=$(sbatch --parsable \
  --dependency="afterok:${pair_job}" \
  --export="ALL,MASLD_PROJECT_ROOT=${PROJECT_ROOT},PLAN45_CANDIDATE_ID=${PLAN45_ORIENTATION_CANDIDATE_ID}" \
  "${SCRIPT_ROOT}/run_orientation_validation.sbatch")
validator_job=${validator_job_raw%%;*}
routing_job_raw=$(sbatch --parsable \
  --dependency="afterok:${validator_job}" \
  --export="ALL,MASLD_PROJECT_ROOT=${PROJECT_ROOT},PLAN45_CANDIDATE_ID=${PLAN45_ROUTING_CANDIDATE_ID},PLAN45_ORIENTATION_ROOT=${orientation_root},PLAN45_LINEAGE_ROOT=${lineage_root},PLAN45_ACCESSIBILITY_ROOT=${accessibility_root}" \
  "${SCRIPT_ROOT}/run_lineage_accessibility_routing.sbatch")
routing_job=${routing_job_raw%%;*}
guide_job_raw=$(sbatch --parsable \
  --dependency="afterok:${routing_job}" \
  --export="ALL,MASLD_PROJECT_ROOT=${PROJECT_ROOT},PLAN45_CANDIDATE_ID=${PLAN45_GUIDE_CANDIDATE_ID},PLAN45_ROUTING_ROOT=${routing_root}" \
  "${SCRIPT_ROOT}/run_guide_design_worklist.sbatch")
guide_job=${guide_job_raw%%;*}

printf 'orientation_tasks=%s\npair_array=%s\norientation_validator=%s\nrouting=%s\nguide_worklist=%s\n' \
  "${n_tasks}" "${pair_job}" "${validator_job}" "${routing_job}" "${guide_job}"
