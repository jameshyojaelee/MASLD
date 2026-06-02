#!/bin/bash
#SBATCH --job-name=tahoe_v2_d2
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=Analysis/Perturbation/scripts/d2_reversal/logs/tahoe_v2_d2_%j.log
#SBATCH --error=Analysis/Perturbation/scripts/d2_reversal/logs/tahoe_v2_d2_%j.log

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

: "${HITS:=Analysis/Perturbation/data/hits/d2_reversal_hits_invocab.csv}"
: "${CONTEXT:=F__Steatohepatitis__Progressor}"
: "${MODALITY:=zero_shot}"
: "${REFERENCE:=ref_a}"
: "${TAHOE_HEPG2_ONLY:=1}"   # v1 pseudobulk-DE: stick to HepG2/C3A shards
: "${TAHOE_ST_MODE:=zeroshot}"
: "${TAHOE_ST_DRY:=0}"

export CUDA_VISIBLE_DEVICES=""
export TAHOE_HEPG2_ONLY
export TAHOE_ST_MODE
export TAHOE_ST_DRY

echo "[$(date)] tahoe_v2_d2 starting: HITS=$HITS CONTEXT=$CONTEXT REFERENCE=$REFERENCE"
echo "  TAHOE_HEPG2_ONLY=$TAHOE_HEPG2_ONLY TAHOE_ST_MODE=$TAHOE_ST_MODE TAHOE_ST_DRY=$TAHOE_ST_DRY"
echo "  SLURM_JOB_ID=$SLURM_JOB_ID  HOST=$(hostname)"

.envs/perturbation_state/bin/python \
  Analysis/Perturbation/scripts/d2_reversal/tahoe_reversal_v2_runner.py \
  --modality "$MODALITY" \
  --context "$CONTEXT" \
  --reference "$REFERENCE" \
  --hits "$HITS"

echo "[$(date)] tahoe_v2_d2 done."
