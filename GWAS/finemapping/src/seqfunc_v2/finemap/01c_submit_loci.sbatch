#!/bin/bash -l
#SBATCH --job-name=SuSiE
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
# Log paths are RELATIVE and resolved against the submission working
# directory, which submit_uniform35_v2.sh sets to $UNIFORM35_RUN_ROOT.
# They were absolute paths into runs/uniform35_v2_2026-07-13/logs/, which
# SLURM evaluates before the script body exports the env var -- so a v3 run
# silently wrote its logs into the v2 run directory.  The submit script also
# passes --output/--error explicitly, which overrides these.
#SBATCH --output=logs/submit_loci_%j.out
#SBATCH --error=logs/submit_loci_%j.err

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
