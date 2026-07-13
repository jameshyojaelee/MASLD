#!/usr/bin/env bash
#SBATCH --job-name=MPRA
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRA_gate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRA_gate_%j.err
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
export MASLD_PROJECT_ROOT=${ROOT}
micromamba run -n rnaseq Rscript "${ROOT}/GWAS/finemapping/src/seqfunc_v2/mpra/04_gate.R"

