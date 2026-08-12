#!/bin/bash
#SBATCH --job-name=proglncrna
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/noncoding_resource/slurm_logs/proglncrna_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/noncoding_resource/slurm_logs/proglncrna_%j.err

set -euo pipefail

: "${IDENTITY_JOB_ID:?Submit with --export=ALL,IDENTITY_JOB_ID=<job_id>}"
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE_ROOT="${PROJECT_ROOT}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11"
IDENTITY="${CANDIDATE_ROOT}/annotation/GENCODE_v49/run_${IDENTITY_JOB_ID}/gene_identity.tsv"
MEMBERSHIP="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
OUTPUT="${CANDIDATE_ROOT}/program_annotation/run_${SLURM_JOB_ID}"

python3 "${PROJECT_ROOT}/scripts/manuscript/noncoding_resource/annotate_program_lncrna_content.py" \
    --membership "${MEMBERSHIP}" \
    --identity "${IDENTITY}" \
    --output "${OUTPUT}"
