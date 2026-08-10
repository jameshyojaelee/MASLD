#!/usr/bin/env bash
#SBATCH --job-name=bridgeget
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --array=0-2
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/bridgeget_%A_%a.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/bridgeget_%A_%a.err

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT="$ROOT/Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge"
SOURCE_GROUPS=(human crop protein)
cd "$ROOT"
python3 "$SCRIPT_ROOT/02_acquire_sources.py" --group "${SOURCE_GROUPS[$SLURM_ARRAY_TASK_ID]}"
