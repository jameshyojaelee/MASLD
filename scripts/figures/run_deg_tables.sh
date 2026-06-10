#!/bin/bash
#SBATCH --job-name=DEGtables
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/deg_tables_%j.log
# cpu is the default; when congested, override to the open io/bigmem fast lane:
#   sbatch --partition=io --qos=interactive --mem=64G run_deg_tables.sh

set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$PROJECT_ROOT"
mkdir -p logs

echo "[$(date)] Starting deg_tables on $(hostname)  JOB=${SLURM_JOB_ID:-local}  CPUs=${SLURM_CPUS_PER_TASK:-4}"
micromamba run -n rnaseq Rscript scripts/figures/figS_disease_signature_sweep_deg_tables.R
echo "[$(date)] Done."
