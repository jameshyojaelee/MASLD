#!/bin/bash -l
# run_finemapping_coloc_targeted.sh
# Orchestrator for the targeted COLOC fine-mapping re-run.
# Submits 02 (per-locus prep) + 03 (SuSiE+CARMA array) chains per GWAS,
# using the lead-SNP files in data/lead_snps_coloc_targeted/.
#
# LD panel: PolyFun UKBB EUR for EUR ancestry; auto-fallback to 1KG for
# EAS/AFR/SAS via finemapping_functions.R (LD_PANEL=polyfun).
#
# Outputs go to output/<study>_coloc_targeted/EUR_0.5Mb/  (parallel namespace
# so the existing fine-mapping run is not clobbered).
#
# Pre-requisite: src/00c_prep_coloc_targeted_leads.R has been run.
#
# Usage:
#   bash run_finemapping_coloc_targeted.sh                # all GWAS with leads
#   bash run_finemapping_coloc_targeted.sh UKBB_GGT       # one GWAS

set -euo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

REGISTRY="config/gwas_registry.tsv"
LEAD_DIR="data/lead_snps_coloc_targeted"

if [ ! -d "${LEAD_DIR}" ]; then
  echo "ERROR: ${LEAD_DIR} does not exist."
  echo "Run: micromamba run -n rnaseq Rscript src/00c_prep_coloc_targeted_leads.R"
  exit 1
fi

STUDY_FILTER="${1:-}"

export LD_PANEL=polyfun
echo "============================================================"
echo "MASLD Fine-Mapping (COLOC-targeted re-run)"
echo "LD_PANEL=${LD_PANEL}  (PolyFun UKBB EUR; 1KG fallback for AFR/EAS/SAS)"
echo "Lead SNPs: ${LEAD_DIR}"
echo "Output namespace: output/<study>_coloc_targeted/"
[ -n "${STUDY_FILTER}" ] && echo "Single study: ${STUDY_FILTER}"
echo "============================================================"

mkdir -p logs

n_submit=0
total_loci=0
tail -n +2 "${REGISTRY}" | while IFS=$'\t' read -r study_name sumstats_path leadsnps_orig ancestry trait_type N_tot N_cases ld_panel window_mb; do
  if [ -n "${STUDY_FILTER}" ] && [ "${study_name}" != "${STUDY_FILTER}" ]; then
    continue
  fi

  TARGETED_LEADS="${LEAD_DIR}/${study_name}_leadSNPs.tsv"
  if [ ! -f "${TARGETED_LEADS}" ]; then
    continue
  fi

  case "${ancestry}" in
    EUR) ld_pop="EUR" ;;
    EAS) ld_pop="EAS" ;;
    AFR) ld_pop="AFR" ;;
    SAS) ld_pop="SAS" ;;
    *) echo "[skip] ${study_name}: unknown ancestry ${ancestry}"; continue ;;
  esac

  TAGGED_STUDY="${study_name}_coloc_targeted"
  N_LOCI=$(( $(wc -l < "${TARGETED_LEADS}") - 1 ))
  if [ "${N_LOCI}" -lt 1 ]; then
    continue
  fi

  echo ""
  echo "--- ${study_name} -> ${TAGGED_STUDY} (${N_LOCI} loci) ---"

  # Stage 1: per-locus sumstats prep (single job)
  PREP_JOB=$(sbatch --parsable \
    --export=ALL,LD_PANEL=${LD_PANEL} \
    --job-name="prep_${TAGGED_STUDY}" \
    --output="logs/prep_${TAGGED_STUDY}_%j.out" \
    --error="logs/prep_${TAGGED_STUDY}_%j.err" \
    src/02_prep_locus_ss.sh \
    "${TAGGED_STUDY}" "${sumstats_path}" "${TARGETED_LEADS}" "${ld_pop}" "${window_mb}")
  echo "  prep job : ${PREP_JOB}"

  # Stage 2+3: SuSiE+CARMA per locus (array, depends on Stage 1)
  FM_JOB=$(sbatch --parsable \
    --export=ALL,LD_PANEL=${LD_PANEL} \
    --dependency=afterok:${PREP_JOB} \
    --array=1-${N_LOCI} \
    --job-name="fm_${TAGGED_STUDY}" \
    --output="logs/fm_${TAGGED_STUDY}_%A_%a.out" \
    --error="logs/fm_${TAGGED_STUDY}_%A_%a.err" \
    src/03_run_fm_per_locus.sh \
    "${TAGGED_STUDY}" "${ld_pop}" "${TARGETED_LEADS}" \
    "${N_tot}" "${N_cases}" "${window_mb}" "${ancestry}")
  echo "  fm  array: ${FM_JOB} (1-${N_LOCI})"

  n_submit=$((n_submit + 1))
  total_loci=$((total_loci + N_LOCI))
done

echo ""
echo "============================================================"
echo "Submitted ${n_submit} GWAS, ${total_loci} fine-mapping array tasks total."
echo "Monitor:  squeue -u \$(whoami) -o '%.10i %.10P %.20j %.8u %.2t %.10M %.6D %R' | grep _coloc_targeted"
echo "============================================================"
