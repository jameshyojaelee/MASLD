#!/usr/bin/env bash
#SBATCH --job-name=MPRA
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=../results/seqfunc/logs/MPRA_%j.out
#SBATCH --error=../results/seqfunc/logs/MPRA_%j.err

set -euo pipefail
ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
mkdir -p "${ROOT}/GWAS/finemapping/results/seqfunc/logs"
bash "${ROOT}/GWAS/finemapping/src/92_download_hu2025_mpra.sh"
export MPRA_NBOOT="${MPRA_NBOOT:-2000}"
micromamba run -n rnaseq Rscript \
  "${ROOT}/GWAS/finemapping/src/92_hu2025_mpra_benchmark.R"
