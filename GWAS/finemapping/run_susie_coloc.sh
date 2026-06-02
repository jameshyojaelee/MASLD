#!/bin/bash -l
# run_susie_coloc.sh — Submit SuSiE-COLOC for all EUR GWAS in the registry
# Reads gwas_registry.tsv, submits 06_susie_coloc.sh array jobs for each EUR study,
# then chains 07_combine_susie_coloc.R after all arrays complete.
#
# Usage:
#   bash run_susie_coloc.sh              # All 24 EUR GWAS
#   bash run_susie_coloc.sh UKBB_ALT     # Single study

set -euo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

REGISTRY="config/gwas_registry.tsv"
STUDY_FILTER="${1:-}"
mkdir -p logs

if [ ! -f "${REGISTRY}" ]; then
    echo "ERROR: Registry file not found: ${REGISTRY}"
    exit 1
fi

echo "============================================================"
echo "SuSiE-COLOC Pipeline — GWAS × Broadaway Liver eQTL"
echo "Registry: ${REGISTRY}"
if [ -n "${STUDY_FILTER}" ]; then
    echo "Running single study: ${STUDY_FILTER}"
fi
echo "============================================================"

COLOC_JOBS=""
N_SUBMITTED=0

# Use process substitution (not pipe) to avoid subshell — variables persist after loop
while IFS=$'\t' read -r study_name sumstats_path leadsnps_path ancestry trait_type N_tot N_cases ld_panel window_mb; do
    # Skip EAS studies (SuSiE-COLOC requires ancestry-matched LD; Broadaway eQTLs are EUR)
    if [ "${ancestry}" != "EUR" ]; then
        echo "  SKIP: ${study_name} (${ancestry} ancestry — EUR eQTLs only)"
        continue
    fi

    # Skip if filtering to a specific study
    if [[ -n "${STUDY_FILTER}" && "${study_name}" != "${STUDY_FILTER}" ]]; then
        continue
    fi

    # Verify summary stats file exists
    if [ ! -f "${sumstats_path}" ]; then
        echo "  WARNING: ${study_name} sumstats not found: ${sumstats_path}, skipping"
        continue
    fi

    echo ""
    echo "--- ${study_name} (${trait_type}, N=${N_tot}) ---"

    JOB_ID=$(sbatch --parsable \
        --job-name="coloc_${study_name}" \
        --output="logs/coloc_${study_name}_%A_%a.out" \
        --error="logs/coloc_${study_name}_%A_%a.err" \
        --export="ALL,GWAS_NAME=${study_name}" \
        src/06_susie_coloc.sh)

    echo "  Array job: ${JOB_ID} (22 tasks, chr1-22)"
    COLOC_JOBS="${COLOC_JOBS}:${JOB_ID}"
    N_SUBMITTED=$((N_SUBMITTED + 1))

done < <(tail -n +2 "${REGISTRY}")

echo ""
echo "============================================================"
echo "Submitted ${N_SUBMITTED} array jobs (${N_SUBMITTED} × 22 = $((N_SUBMITTED * 22)) tasks)"

# Submit combiner as dependency after all array jobs complete
if [ -n "${COLOC_JOBS}" ]; then
    COMBINE_JOB=$(sbatch --parsable \
        --dependency=afterany${COLOC_JOBS} \
        --partition=cpu \
        --cpus-per-task=1 \
        --mem=32G \
        --time=48:00:00 \
        --job-name="coloc_combine" \
        --output="logs/coloc_combine_%j.out" \
        --error="logs/coloc_combine_%j.err" \
        --wrap="bash -c 'cd ${FM_DIR} && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate finemapping && Rscript src/07_combine_susie_coloc.R'")

    echo "Combiner job: ${COMBINE_JOB} (runs after all arrays)"
fi
echo "Monitor: squeue -u \$(whoami) | grep coloc"
echo "============================================================"
