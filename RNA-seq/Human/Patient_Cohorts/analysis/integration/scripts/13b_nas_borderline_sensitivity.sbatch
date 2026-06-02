#!/bin/bash
#SBATCH --job-name=dream-NAS
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=logs/13b_nas_borderline_%j.out
#SBATCH --error=logs/13b_nas_borderline_%j.err

set -eo pipefail

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

echo "=== Job started: $(date) ==="
echo "Node: ${SLURMD_NODENAME:-unknown}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE:-unknown} MB"

Rscript "${SCRIPTS}/13b_nas_borderline_sensitivity.R"

echo "=== Job finished: $(date) ==="
