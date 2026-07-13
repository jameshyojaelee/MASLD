#!/bin/bash -l
#SBATCH --job-name=Rscript
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=GWAS/MR_Data/MVP/logs/replication_lookup_%j.out
#SBATCH --error=GWAS/MR_Data/MVP/logs/replication_lookup_%j.err

# Cross-ancestry variant-level replication lookup (EUR-discovered leads -> non-EUR MVP arms).
# Reads full non-EUR sumstats (up to ~2.2 GB) into memory -> 128 GB is generous.
#
# Usage:
#   sbatch scripts/mvp/run_replication_lookup.sh
#   REPL_PHENOS=ALL sbatch scripts/mvp/run_replication_lookup.sh
#   REPL_PROXY_KB=0 sbatch scripts/mvp/run_replication_lookup.sh   # exact-only

set -eo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "${BASE}"

eval "$(micromamba shell hook -s bash)"
set +u                          # micromamba activate references unset vars
micromamba activate rnaseq
set -u

echo "=== Cross-ancestry replication lookup ==="
echo "Host: $(hostname)  Start: $(date)"
echo "PHENOS=${REPL_PHENOS:-<default>}  PROXY_KB=${REPL_PROXY_KB:-<default>}  CLUMP_KB=${REPL_CLUMP_KB:-<default>}"

Rscript GWAS/MR_Data/MVP/cross_ancestry_replication_lookup.R

echo "=== Done: $(date) ==="
