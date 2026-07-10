#!/bin/bash
#SBATCH --job-name=fgsea
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/14g_grouped_nas_gsea_%j.out
#SBATCH --error=logs/14g_grouped_nas_gsea_%j.err

# drop `-u`: micromamba binutils activate.d references $ADDR2LINE unbound
set -eo pipefail

SCRIPTS="${SLURM_SUBMIT_DIR:-$(dirname "$(realpath "$0")")}"
BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
mkdir -p "${SCRIPTS}/logs"

echo "=== Job started: $(date) ==="
echo "Node: ${SLURMD_NODENAME:-unknown}   BASE=${BASE}"

# 1) fgsea on the grouped-NAS contrasts -> nas_grouped_gsea.csv
micromamba run -n rnaseq Rscript "${SCRIPTS}/14g_grouped_nas_gsea.R"

# 2) render the NAS-activity Hallmark escalation dose-response figure
micromamba run -n rnaseq Rscript "${BASE}/scripts/figures/nas_activity_hallmark_escalation.R"

echo "=== Job finished: $(date) ==="
