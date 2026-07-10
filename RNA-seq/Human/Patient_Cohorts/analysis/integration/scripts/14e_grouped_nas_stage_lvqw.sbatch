#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=logs/14e_grouped_nas_stage_lvqw_%j.out
#SBATCH --error=logs/14e_grouped_nas_stage_lvqw_%j.err

# NOTE: drop `-u` — micromamba's binutils activate.d references $ADDR2LINE
# unbound, which trips `set -u` before the R script ever runs.
set -eo pipefail

# Use the submission dir under SLURM (realpath "$0" points at the unwritable
# node spool copy on some partitions, e.g. io); fall back to $0 for local runs.
SCRIPTS="${SLURM_SUBMIT_DIR:-$(dirname "$(realpath "$0")")}"
mkdir -p "${SCRIPTS}/logs"

echo "=== Job started: $(date) ==="
echo "Node: ${SLURMD_NODENAME:-unknown}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

micromamba run -n rnaseq Rscript "${SCRIPTS}/14e_grouped_nas_stage_lvqw.R"

echo "=== Job finished: $(date) ==="
