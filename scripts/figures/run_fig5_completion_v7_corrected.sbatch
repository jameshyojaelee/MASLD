#!/bin/bash
#SBATCH --job-name=fig5_v7_corr
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig5_v7_corr_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig5_v7_corr_%j.err

# Figure 5 v7 candidate on the corrected inputs (review items 2, 7, 8):
#   5C from the corrected-control F_five fit (BULK_DEG_FILE),
#   5D from the COLOC-eligible primary pairs written by run_fig5d_coloc_r2.sbatch.
# Reads the 5D census at run time, then runs run_fig5_completion_candidate.sbatch
# unchanged in this allocation.
#   CANDIDATE_ROOT  new figures/candidates/fig5-completion-<date>-v7 (must not exist)
#   BULK_DEG_FILE   corrected deg_results.csv
#   FIG5D_CONTEXT   <atac-context-v3-coloc-r2-*>/genetics/context
set -euo pipefail
: "${CANDIDATE_ROOT:?}"; : "${BULK_DEG_FILE:?}"; : "${FIG5D_CONTEXT:?}"
test ! -e "$CANDIDATE_ROOT"
test -s "$BULK_DEG_FILE"
census() { awk -F'\t' -v m="$1" '$1 == m {print $3}' "$FIG5D_CONTEXT/restriction_summary.tsv"; }
export FIG5D_PRIMARY_PAIRS="$FIG5D_CONTEXT/genetic_lineage_context_primary_pairs.tsv"
export FIG5D_EXPECTED_ROWS=$(census lineage_rows)
export FIG5D_EXPECTED_PAIRS=$(census gene_study_pairs)
export FIG5D_EXPECTED_GENES=$(census genes)
test -s "$FIG5D_PRIMARY_PAIRS"
test -n "$FIG5D_EXPECTED_ROWS" && test -n "$FIG5D_EXPECTED_PAIRS" && test -n "$FIG5D_EXPECTED_GENES"
echo "FIG5D census: rows=$FIG5D_EXPECTED_ROWS pairs=$FIG5D_EXPECTED_PAIRS genes=$FIG5D_EXPECTED_GENES"
export CANDIDATE_ROOT BULK_DEG_FILE
export FIG5B_BULK_DEG="$BULK_DEG_FILE"   # render 5B from the corrected fit instead of copying it
exec bash scripts/figures/run_fig5_completion_candidate.sbatch
