#!/bin/bash
#SBATCH --job-name=studyaf
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/study_af/merge_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/study_af/merge_%j.err

# Collects the per-study fragments into config/gwas_af_sources.tsv.  Exits
# non-zero if any of the 35 studies is still missing a fragment.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
Rscript "${PROJECT_ROOT}/GWAS/finemapping/src/seqfunc_v2/finemap/00b_build_study_af.R" --merge
