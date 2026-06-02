#!/bin/bash -l
# run_finemapping.sh — Master orchestrator for MASLD fine-mapping pipeline
# Reads study list from config/gwas_registry.tsv and runs:
#   Stage 1: Prep per-locus summary stats
#   Stage 2+3: Fine-mapping (SuSiE + CARMA) per locus via SLURM array
#
# Usage: bash run_finemapping.sh [study_name]
#   If study_name is provided, only run that study
#   Otherwise, run all studies in the registry

set -euo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

REGISTRY="config/gwas_registry.tsv"

if [ ! -f "${REGISTRY}" ]; then
    echo "ERROR: Registry file not found: ${REGISTRY}"
    exit 1
fi

# Filter to specific study if provided
STUDY_FILTER="${1:-}"

echo "============================================================"
echo "MASLD Fine-Mapping Pipeline"
echo "Registry: ${REGISTRY}"
if [ -n "${STUDY_FILTER}" ]; then
    echo "Running single study: ${STUDY_FILTER}"
fi
echo "============================================================"

# Read registry (TSV: study_name, sumstats_path, leadsnps_path, ancestry, trait_type, N_tot, N_cases, ld_panel, window_mb)
tail -n +2 "${REGISTRY}" | while IFS=$'\t' read -r study_name sumstats_path leadsnps_path ancestry trait_type N_tot N_cases ld_panel window_mb; do
    # Skip if filtering to a specific study
    if [ -n "${STUDY_FILTER}" ] && [ "${study_name}" != "${STUDY_FILTER}" ]; then
        continue
    fi

    # Determine LD population label
    if [ "${ancestry}" == "EUR" ]; then
        ld_pop="EUR"
    elif [ "${ancestry}" == "EAS" ]; then
        ld_pop="EAS"
    else
        echo "WARNING: Unknown ancestry ${ancestry} for ${study_name}, skipping"
        continue
    fi

    # Check if EAS LD panel exists
    if [ "${ancestry}" == "EAS" ] && [ ! -d "data/ld_ref/1kg_eas/chr1" ]; then
        echo "WARNING: EAS LD panel not built yet, skipping ${study_name}"
        continue
    fi

    # Count lead SNPs (subtract header)
    N_LOCI=$(( $(wc -l < "${leadsnps_path}") - 1 ))
    if [ "${N_LOCI}" -lt 1 ]; then
        echo "WARNING: No lead SNPs found for ${study_name}, skipping"
        continue
    fi

    echo ""
    echo "--- Study: ${study_name} ---"
    echo "  Ancestry: ${ancestry}, N_tot: ${N_tot}, N_cases: ${N_cases}"
    echo "  Lead SNPs: ${N_LOCI}, Window: ${window_mb} Mb"
    echo "  LD panel: ${ld_panel}"

    # Stage 1: Prep per-locus summary stats
    echo "  Submitting Stage 1 (prep locus SS)..."
    PREP_JOB=$(sbatch --parsable \
        --job-name="prep_${study_name}" \
        --output="logs/prep_${study_name}_%j.out" \
        --error="logs/prep_${study_name}_%j.err" \
        src/02_prep_locus_ss.sh \
        "${study_name}" "${sumstats_path}" "${leadsnps_path}" "${ld_pop}" "${window_mb}")
    echo "  Stage 1 job: ${PREP_JOB}"

    # Stage 2+3: Fine-mapping per locus (array job, depends on Stage 1)
    echo "  Submitting Stage 2+3 (fine-mapping, ${N_LOCI} loci)..."
    FM_JOB=$(sbatch --parsable \
        --dependency=afterok:${PREP_JOB} \
        --array=1-${N_LOCI} \
        --job-name="fm_${study_name}" \
        src/03_run_fm_per_locus.sh \
        "${study_name}" "${ld_pop}" "${leadsnps_path}" \
        "${N_tot}" "${N_cases}" "${window_mb}" "${ancestry}")
    echo "  Stage 2+3 job array: ${FM_JOB}"

done

echo ""
echo "============================================================"
echo "All jobs submitted. Monitor with: squeue -u \$(whoami) | grep fm_"
echo "============================================================"
