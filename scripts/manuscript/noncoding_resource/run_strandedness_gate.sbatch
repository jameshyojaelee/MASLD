#!/bin/bash
#SBATCH --job-name=strandgate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/noncoding_resource/slurm_logs/strandgate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/noncoding_resource/slurm_logs/strandgate_%j.err

set -euo pipefail

: "${STRAND_RUN_JOB_ID:?Submit with --export=ALL,STRAND_RUN_JOB_ID=<job_id>}"
: "${STRAND_POLICY:=all_samples_fourfold_v1}"
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE_ROOT="${PROJECT_ROOT}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11"
FEATURECOUNTS_RUN="${CANDIDATE_ROOT}/strand_audit/GSE162694/run_${STRAND_RUN_JOB_ID}"
OUTPUT_DIR="${CANDIDATE_ROOT}/source_gate/run_${SLURM_JOB_ID}"
METADATA="${PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"

if [[ ! -f "${FEATURECOUNTS_RUN}/COMPLETE" ]]; then
    echo "Upstream featureCounts audit is not complete: ${FEATURECOUNTS_RUN}" >&2
    exit 1
fi

python3 "${PROJECT_ROOT}/scripts/manuscript/noncoding_resource/audit_strandedness.py" \
    --project-root "${PROJECT_ROOT}" \
    --featurecounts-run "${FEATURECOUNTS_RUN}" \
    --policy "${STRAND_POLICY}" \
    --metadata "${METADATA}" \
    --output-dir "${OUTPUT_DIR}"

python3 -m pip freeze > "${OUTPUT_DIR}/pip-freeze.txt"
sha256sum "${OUTPUT_DIR}"/*.tsv > "${OUTPUT_DIR}/outputs.sha256"
touch "${OUTPUT_DIR}/COMPLETE"
