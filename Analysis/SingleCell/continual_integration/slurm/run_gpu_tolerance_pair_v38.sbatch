#!/usr/bin/env bash
#SBATCH --job-name=masld-gpu-tol
#SBATCH --partition=gpu
#SBATCH --qos=interactive
#SBATCH --gres=gpu:b6k:2
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/gpu-tolerance-v38-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/gpu-tolerance-v38-%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
SPECS="${PIPELINE_ROOT}/candidates/gpu-tolerance-specs-v38-2026-08-15T163000Z"

bash "${PIPELINE_ROOT}/candidates/interactive-v100-pair-runner-v1.sh" \
  "${SPECS}/replica_a.json" "${SPECS}/replica_b.json"
