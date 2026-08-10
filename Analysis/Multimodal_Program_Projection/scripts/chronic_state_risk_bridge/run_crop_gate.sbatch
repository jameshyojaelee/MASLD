#!/usr/bin/env bash
#SBATCH --job-name=cropgate
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/cropgate_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/cropgate_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$ROOT"
python3 Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/05_gate_cropseq.py
