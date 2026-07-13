#!/bin/bash
#SBATCH --job-name=noneurAudit
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/noneur_audit_%j.log

set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$MASLD_PROJECT_ROOT"
echo "host=$(hostname) start=$(date)"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/figures/audit_noneur_gws.R
echo "end=$(date)"
