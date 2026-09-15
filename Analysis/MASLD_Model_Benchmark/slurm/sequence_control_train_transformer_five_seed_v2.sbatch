#!/usr/bin/env bash
#SBATCH --job-name=model-training-704
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
: "${OUTER_FOLD:=0}"
: "${LINEAGE:=hepatocyte}"
: "${INPUT_JOB_ID:=21064244}"
: "${INPUT_ARTIFACTS_SHA256:=fc26ec558ea51f74212bea120dab59053e5f2856494407bb0a717e9eaecc0208}"
SPLIT_ID="donor${OUTER_FOLD}_genomic${OUTER_FOLD}"
INPUT_PATH="${ROOT}/executions/chrombpnet-${LINEAGE}-${SPLIT_ID}-inputs-${INPUT_JOB_ID}"
if [[ "${LINEAGE}:${OUTER_FOLD}:${INPUT_JOB_ID}" == hepatocyte:0:21064244 ]]; then
  INPUT_PATH="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-inputs-21064244"
fi
WRAPPER="${ROOT}/slurm/sequence_control_train_transformer_five_seed_v2.sbatch"
exec bash "${ROOT}/scripts/sequence_control_production_driver_five_seed_v2.sh" \
  sequence_transformer_control \
  "${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993" \
  3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b \
  "${INPUT_PATH}" \
  "${INPUT_ARTIFACTS_SHA256}" \
  32 \
  "${WRAPPER}"
