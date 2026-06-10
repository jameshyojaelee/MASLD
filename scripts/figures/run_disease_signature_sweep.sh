#!/bin/bash
#SBATCH --job-name=DEGsweep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/disease_signature_sweep_%j.log
# NOTE: cpu is the default; when it is congested, submit with an override to the
# open io/bigmem fast lane (bypasses the nslab 7TB cap, max 4 interactive jobs):
#   sbatch --partition=io --qos=interactive --mem=64G run_disease_signature_sweep.sh

set -euo pipefail

# Under SLURM, $0 is a spool copy — resolve the repo root from the known absolute
# path (same default the R script uses) so cd/mkdir land in a writable place.
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$PROJECT_ROOT"
mkdir -p logs

echo "[$(date)] Starting disease_signature_sweep on $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-local}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}, MEM: 64G"

micromamba run -n rnaseq Rscript scripts/figures/figS_disease_signature_sweep.R

echo "[$(date)] Done."
