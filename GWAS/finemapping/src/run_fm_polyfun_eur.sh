#!/bin/bash
# run_fm_polyfun_eur.sh — Phase 9e PolyFun standalone fine-mapping
# (SuSiE + CARMA per locus across all 17 EUR GWAS).
#
# Submitter (login-node script). Loops over EUR GWAS and submits one
# sbatch array per GWAS to the cluster. ~170 loci total.
#
# Run:
#   bash run_fm_polyfun_eur.sh
#
# Per-locus task:
#   - Job name: polyfun-fm
#   - 16 GB RAM, 4 CPUs, 48h wall
#   - Partitions: cpu, io, gpu_interactive (multi-partition spread)
#   - CARMA-bound: 4-6 hr typical
#   - Idempotent skip: 03_run_fm_per_locus.R checks for existing output
#
# Env activation handled inside run_fm_polyfun_eur_jobwrapper.sh
# (uses absolute paths to micromamba, so portable across users).
#
# Outputs: output/<gwas>/polyfun_eur_0.5Mb/{susie,CARMA}/

set -o pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

NEW_LD_POP="polyfun_eur"
REGISTRY="config/gwas_registry.tsv"
JOBNAME="polyfun-fm"
PARTITIONS="cpu,io,gpu_interactive"
JOBWRAPPER="src/run_fm_polyfun_eur_jobwrapper.sh"

echo "============================================================"
echo "Phase 9e — PolyFun FM (SuSiE+CARMA), all 17 EUR GWAS"
echo "Job name:   ${JOBNAME}"
echo "Partitions: ${PARTITIONS}"
echo "============================================================"

SUBMITTED=0
SKIPPED=0

tail -n +2 "${REGISTRY}" \
  | while IFS=$'\t' read -r study_name sumstats_path leadsnps_path ancestry trait_type N_tot N_cases ld_panel window_mb; do
    [ "${ancestry}" != "EUR" ] && continue

    win="${window_mb:-0.5}"
    OLD_SS_DIR="output/${study_name}/EUR_${win}Mb/ss"
    NEW_BASE="output/${study_name}/${NEW_LD_POP}_${win}Mb"
    NEW_SS_DIR="${NEW_BASE}/ss"
    mkdir -p "${NEW_SS_DIR}"

    n_linked=0
    if [ -d "${OLD_SS_DIR}" ]; then
      for f in "${OLD_SS_DIR}"/*.txt; do
        [ -f "$f" ] || continue
        dest="${NEW_SS_DIR}/$(basename "$f")"
        [ -e "$dest" ] && continue
        ln -s "$(realpath "$f")" "$dest" 2>/dev/null && n_linked=$((n_linked+1))
      done
    else
      echo "  [skip] ${study_name}: no per-locus ss/ at ${OLD_SS_DIR}"
      SKIPPED=$((SKIPPED+1))
      continue
    fi

    N_LOCI=$(( $(wc -l < "${leadsnps_path}") - 1 ))
    if [ "${N_LOCI}" -lt 1 ]; then
      echo "  [skip] ${study_name}: 0 lead SNPs"
      SKIPPED=$((SKIPPED+1))
      continue
    fi

    echo ""
    echo "--- ${study_name}: ${N_LOCI} loci, ${n_linked} ss/ links ---"

    JOB=$(sbatch --parsable \
      --export=ALL,LD_PANEL=polyfun \
      --array=1-${N_LOCI} \
      --job-name="${JOBNAME}" \
      --partition=${PARTITIONS} \
      --cpus-per-task=4 \
      --mem=64G \
      --time=90:00:00 \
      --output="logs/polyfun_fm_${study_name}_%A_%a.out" \
      --error="logs/polyfun_fm_${study_name}_%A_%a.err" \
      ${JOBWRAPPER} \
      "${study_name}" "${NEW_LD_POP}" "${leadsnps_path}" \
      "${N_tot}" "${N_cases}" "${win}" "EUR")
    echo "  submitted: ${JOB}  (array 1-${N_LOCI})"
    SUBMITTED=$((SUBMITTED+1))
  done

echo ""
echo "============================================================"
echo "Submitted ${SUBMITTED} per-GWAS array jobs (${SKIPPED} skipped)."
echo "Monitor: squeue -u \$(whoami) -n ${JOBNAME}"
echo "Outputs: output/{study}/${NEW_LD_POP}_0.5Mb/{susie,CARMA}/"
echo "============================================================"
