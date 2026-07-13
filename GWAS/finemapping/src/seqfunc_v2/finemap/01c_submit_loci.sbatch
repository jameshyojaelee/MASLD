#!/bin/bash -l
#SBATCH --job-name=SuSiE
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/runs/uniform35_v2_2026-07-13/logs/submit_loci_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/runs/uniform35_v2_2026-07-13/logs/submit_loci_%j.err

set -euo pipefail
ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
export MASLD_PROJECT_ROOT="${ROOT}"
export UNIFORM35_RUN_ROOT="${UNIFORM35_RUN_ROOT:-${ROOT}/GWAS/finemapping/runs/uniform35_v2_2026-07-13}"
SRC="${ROOT}/GWAS/finemapping/src/seqfunc_v2/finemap"

# Finalization is the hard preparation gate. Refuse to submit an empty or
# partially materialized locus portfolio.
test -s "${UNIFORM35_RUN_ROOT}/config/locus_manifest.tsv"
test -s "${UNIFORM35_RUN_ROOT}/config/preparation_summary.json"
exec bash "${SRC}/submit_uniform35_v2.sh" run
