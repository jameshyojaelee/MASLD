#!/bin/bash
#SBATCH --job-name=tahoe_v2_d1
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=Analysis/Perturbation/scripts/d1_mechanism/logs/tahoe_v2_d1_%j.log
#SBATCH --error=Analysis/Perturbation/scripts/d1_mechanism/logs/tahoe_v2_d1_%j.log

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Defaults — override via env vars / sbatch --export
: "${HITS:=Analysis/Perturbation/data/hits/d1_mechanism_hits_invocab.csv}"
: "${CONTEXT:=all}"
: "${MODALITY:=zero_shot}"
: "${TAHOE_N_SHARDS:=50}"   # v1 pseudobulk-DE shard scan cap
: "${TAHOE_ST_MODE:=zeroshot}"
: "${TAHOE_ST_DRY:=0}"

export CUDA_VISIBLE_DEVICES=""
export TAHOE_N_SHARDS
export TAHOE_ST_MODE
export TAHOE_ST_DRY

echo "[$(date)] tahoe_v2_d1 starting: HITS=$HITS CONTEXT=$CONTEXT MODALITY=$MODALITY"
echo "  TAHOE_N_SHARDS=$TAHOE_N_SHARDS TAHOE_ST_MODE=$TAHOE_ST_MODE TAHOE_ST_DRY=$TAHOE_ST_DRY"
echo "  SLURM_JOB_ID=$SLURM_JOB_ID  HOST=$(hostname)"

.envs/perturbation_state/bin/python \
  Analysis/Perturbation/scripts/d1_mechanism/tahoe_v2_runner.py \
  --modality "$MODALITY" \
  --context "$CONTEXT" \
  --hits "$HITS"

echo "[$(date)] tahoe_v2_d1 done."
