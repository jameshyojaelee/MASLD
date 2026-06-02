#!/bin/bash -l
# run_fm_topld_nonEUR.sh
# Re-run standalone GWAS fine-mapping (SuSiE + CARMA) for the 11 non-EUR GWAS
# (5 EAS + 3 AFR + 3 SAS) with LD_PANEL=topld so each ancestry's TOP-LD panel
# is used (topld_eas / topld_afr / topld_sas). The existing v1 1KG-LD outputs at
# output/{gwas}/{ANCESTRY}_0.5Mb/ are the baseline for comparison.
#
# Strategy mirrors run_fm_topld_eur.sh:
#   - Reuse existing v1 per-locus summary stat files under
#     output/{gwas}/{ANCESTRY}_{win}Mb/ss/  (prep is LD-independent).
#   - Symlink ss/*.txt into output/{gwas}/topld_<lower(ancestry)>_{win}Mb/ss/.
#   - Submit per-study array via 03_run_fm_per_locus.sh with
#     LD_PANEL=topld exported. EAS_LD_DIR / AFR_LD_DIR / SAS_LD_DIR are NOT
#     pinned, so finemapping_functions.R::get_ld_base_dir() routes the ancestry
#     to its topld_<ancestry>/ panel.
#   - 03_run_fm_per_locus.R needs only .ld + .bim per block (no .bed/.fam),
#     which the topld_<ancestry>/ panels already have — no backfill required.

set -o pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

REGISTRY="config/gwas_registry.tsv"
WIN_MB=0.5
STUDY_FILTER="${1:-}"

echo "============================================================"
echo "MASLD Fine-Mapping — TOP-LD non-EUR (EAS/AFR/SAS)"
echo "Registry: ${REGISTRY}"
echo "Target: 11 non-EUR GWAS, LD_PANEL=topld, output dir suffix=topld_<ancestry>"
if [ -n "${STUDY_FILTER}" ]; then
    echo "Running single study: ${STUDY_FILTER}"
fi
echo "============================================================"

export LD_PANEL=topld
# Explicitly leave EAS_LD_DIR/AFR_LD_DIR/SAS_LD_DIR UNSET so dispatch resolves
# to topld_<ancestry>/ via the panel switch in get_ld_base_dir().
unset EAS_LD_DIR AFR_LD_DIR SAS_LD_DIR

tail -n +2 "${REGISTRY}" | while IFS=$'\t' read -r study_name sumstats_path leadsnps_path ancestry trait_type N_tot N_cases ld_panel window_mb; do
    case "${ancestry}" in EAS|AFR|SAS) ;; *) continue;; esac
    [ -n "${STUDY_FILTER}" ] && [ "${study_name}" != "${STUDY_FILTER}" ] && continue
    win="${window_mb:-0.5}"

    # Lowercase ancestry for output dir
    anc_lower=$(echo "${ancestry}" | tr '[:upper:]' '[:lower:]')
    NEW_LD_POP="topld_${anc_lower}"

    # Existing ss dir from the v1 1KG-LD non-EUR run (e.g., output/<gwas>/EAS_0.5Mb/ss)
    OLD_SS_DIR="output/${study_name}/${ancestry}_${win}Mb/ss"
    NEW_BASE="output/${study_name}/${NEW_LD_POP}_${win}Mb"
    NEW_SS_DIR="${NEW_BASE}/ss"
    mkdir -p "${NEW_SS_DIR}"

    # Symlink each ss .txt if not already linked
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
        echo "WARNING: ${study_name} has 0 lead SNPs (${ancestry}); skipping"
        continue
    fi

    echo ""
    echo "--- ${study_name} (${ancestry}, ${N_LOCI} loci, ${n_linked} ss/ links) ---"

    JOB=$(sbatch --parsable \
        --export=ALL,LD_PANEL=topld \
        --array=1-${N_LOCI} \
        --job-name="fmtopld_${study_name}" \
        --partition=cpu \
        --cpus-per-task=4 \
        --mem=64G \
        --time=48:00:00 \
        --output="logs/fmtopld_${study_name}_%A_%a.out" \
        --error="logs/fmtopld_${study_name}_%A_%a.err" \
        src/03_run_fm_per_locus.sh \
        "${study_name}" "${NEW_LD_POP}" "${leadsnps_path}" \
        "${N_tot}" "${N_cases}" "${win}" "${ancestry}")
    echo "  Submitted: ${JOB}"
done

echo ""
echo "============================================================"
echo "All non-EUR fine-mapping jobs submitted. Monitor: squeue -u \$(whoami) | grep fmtopld_"
echo "Output structure: output/{study}/topld_<ancestry>_*/  (parallel to {ANCESTRY}_0.5Mb/)"
echo "============================================================"
