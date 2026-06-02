#!/bin/bash
#SBATCH --job-name=14b_vs_ctrl
#SBATCH --partition=bigmem
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=logs/14b_stage_vs_healthy_dream_%j.out
#SBATCH --error=logs/14b_stage_vs_healthy_dream_%j.err

set -euo pipefail

SCRIPTS="$(dirname "$(realpath "$0")")"
mkdir -p "${SCRIPTS}/logs"

echo "=== Job started: $(date) ==="
echo "Node: ${SLURMD_NODENAME:-unknown}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

micromamba run -n rnaseq Rscript "${SCRIPTS}/14b_stage_vs_healthy_dream.R"

echo "=== Job finished: $(date) ==="
