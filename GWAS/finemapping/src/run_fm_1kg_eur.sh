#!/bin/bash -l
# run_fm_1kg_eur.sh
# Re-run standalone GWAS fine-mapping (SuSiE + CARMA) for all 17 EUR GWAS
# using the new 1KG EUR LD panel instead of sghatan UKBB EUR.
#
# Strategy:
#   - Reuse existing per-locus summary stat files under
#     output/{gwas}/EUR_{win}Mb/ss/ (prep is LD-independent).
#   - Submit 03_run_fm_per_locus.sh per study with ld_pop="1kg_eur" AND
#     LD_PANEL=1kg exported, so finemapping_functions.R::get_ld_base_dir()
#     resolves to data/ld_ref/1kg_eur/.
#   - Output lands at output/{gwas}/1kg_eur_{win}Mb/ (new dir, does not
#     overwrite the v1 UKBB-based output at output/{gwas}/EUR_{win}Mb/).
# The SuSiE-COLOC production re-run (Phase 4) is a separate pipeline (06_susie_coloc.R).

set -o pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

REGISTRY="config/gwas_registry.tsv"
WIN_MB=0.5
NEW_LD_POP="1kg_eur"

STUDY_FILTER="${1:-}"

echo "============================================================"
echo "MASLD Fine-Mapping — 1KG EUR re-run (sghatan severance)"
echo "Registry: ${REGISTRY}"
echo "Target: all EUR GWAS, LD_PANEL=1kg, output dir suffix=${NEW_LD_POP}"
if [ -n "${STUDY_FILTER}" ]; then
    echo "Running single study: ${STUDY_FILTER}"
fi
echo "============================================================"

# Ensure env var propagates to sbatch child jobs
export LD_PANEL=1kg

# Parse registry (TSV header + 17 EUR rows we care about)
tail -n +2 "${REGISTRY}" | while IFS=$'\t' read -r study_name sumstats_path leadsnps_path ancestry trait_type N_tot N_cases ld_panel window_mb; do
    [ "${ancestry}" != "EUR" ] && continue
    [ -n "${STUDY_FILTER}" ] && [ "${study_name}" != "${STUDY_FILTER}" ] && continue
    # window_mb from registry (fall back to 0.5)
    win="${window_mb:-0.5}"

    # Existing ss dir from the v1 UKBB run (LD-independent; reuse to avoid re-prep)
    OLD_SS_DIR="output/${study_name}/EUR_${win}Mb/ss"
    NEW_BASE="output/${study_name}/${NEW_LD_POP}_${win}Mb"
    NEW_SS_DIR="${NEW_BASE}/ss"
    mkdir -p "${NEW_SS_DIR}"

    # Symlink each .txt if not already linked
    n_linked=0
    if [ -d "${OLD_SS_DIR}" ]; then
        for f in "${OLD_SS_DIR}"/*.txt; do
            [ -f "$f" ] || continue
            dest="${NEW_SS_DIR}/$(basename "$f")"
            [ -e "$dest" ] && continue
            ln -s "$(realpath "$f")" "$dest" 2>/dev/null && n_linked=$((n_linked+1))
        done
    else
        echo "  WARNING: ${OLD_SS_DIR} missing; cannot reuse per-locus SS for ${study_name}"
        continue
    fi

    N_LOCI=$(( $(wc -l < "${leadsnps_path}") - 1 ))
    if [ "${N_LOCI}" -lt 1 ]; then
        echo "WARNING: No lead SNPs for ${study_name}, skipping"
        continue
    fi

    echo ""
    echo "--- ${study_name} (EUR, ${N_LOCI} loci, ${n_linked} ss/ files linked) ---"

    # Submit Stage 2+3 array job with LD_PANEL=1kg, writing to 1kg_eur_Xmb dir
    # Using sbatch --export=ALL,LD_PANEL=1kg so the env var propagates
    JOB=$(sbatch --parsable \
        --export=ALL,LD_PANEL=1kg \
        --array=1-${N_LOCI} \
        --job-name="fm1kg_${study_name}" \
        --partition=cpu \
        --cpus-per-task=4 \
        --mem=64G \
        --time=48:00:00 \
        --output="logs/fm1kg_${study_name}_%A_%a.out" \
        --error="logs/fm1kg_${study_name}_%A_%a.err" \
        src/03_run_fm_per_locus.sh \
        "${study_name}" "${NEW_LD_POP}" "${leadsnps_path}" \
        "${N_tot}" "${N_cases}" "${win}" "EUR")
    echo "  Submitted: ${JOB}"
done

echo ""
echo "============================================================"
echo "All EUR fine-mapping jobs submitted. Monitor with: squeue -u \$(whoami) | grep fm1kg_"
echo "Output structure: output/{study}/1kg_eur_*/  (parallel to EUR_0.5Mb/)"
echo "============================================================"
