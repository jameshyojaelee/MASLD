#!/usr/bin/env bash
#SBATCH --job-name=bridgegate
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/bridgegate_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/bridgegate_%j.err

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS="$ROOT/Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge"
cd "$ROOT"
python3 "$SCRIPTS/03_gate_human_pairs.py"
micromamba run -n rnaseq Rscript "$SCRIPTS/04_fingerprint_human_overlap.R"
