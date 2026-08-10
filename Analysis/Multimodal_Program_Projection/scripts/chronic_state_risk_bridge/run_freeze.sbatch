#!/usr/bin/env bash
#SBATCH --job-name=bridgefreeze
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/bridgefreeze_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/bridgefreeze_%j.err

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT="$ROOT/Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge"
cd "$ROOT"
python3 "$SCRIPT_ROOT/00_freeze_and_correct.py"
python3 "$SCRIPT_ROOT/01_validate_seal.py"
