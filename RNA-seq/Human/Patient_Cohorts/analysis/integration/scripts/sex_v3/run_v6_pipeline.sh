#!/bin/bash
# ---------------------------------------------------------------------------
# sex_v6 reproducibility wrapper
# ---------------------------------------------------------------------------
# Usage:
#   bash run_v6_pipeline.sh --from-step N [--dry-run] [--n-perm N]
#                           [--n-boot N] [--n-seeds N] [--gap-fill]
#
# Steps (per plan §Dependency Graph):
#   1. P1 random-slope refit
#   2. P3 selection-bias sim       \ independent of P1
#   3. P5 50-seed calibration       \ submit at t=0
#   4. P6 Cochran Q                 /
#   5. P2 perm-FDR        (afterok:P1)
#   6. P4 dream-bootstrap (afterok:P1)
#   7. P7 threshold sweep (afterok:P1)
#   8. 8A male char + 8B LOCO magnitude (afterany:P6 for Q col)
#   9. P9 manifest emit + consensus 18 (afterany ALL)
#
# Validates input file existence, captures all JIDs, writes
# pipeline_manifest_v6.json (initial scaffold; final mutation by
# v6_manifest_emit.R).
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"

# Contrast routing (added 2026-05-18 for stage-stratified extension).
# CONTRAST_NAME ∈ {disease_vs_ctrl (default), mash_vs_masl, masl_vs_ctrl, mash_vs_ctrl}.
# Default preserves the legacy `sex_v3/intermediates/` layout; non-default
# contrasts use a sibling subdir `sex_v3/contrast_<name>/intermediates/`.
CONTRAST_NAME="${CONTRAST_NAME:-disease_vs_ctrl}"
if [[ "$CONTRAST_NAME" == "disease_vs_ctrl" ]]; then
  SEXV3_OUT="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3"
else
  SEXV3_OUT="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/contrast_${CONTRAST_NAME}"
fi
INTERMED_DIR="${SEXV3_OUT}/intermediates"
INTEGR_DIR="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
MANIFEST="${SEXV3_OUT}/pipeline_manifest_v6.json"
mkdir -p "${INTERMED_DIR}" "${SEXV3_OUT}/logs"
echo "Pipeline driver: CONTRAST_NAME=${CONTRAST_NAME}"
echo "  SEXV3_OUT     = ${SEXV3_OUT}"
echo "  INTERMED_DIR  = ${INTERMED_DIR}"
echo

# Export so every sbatch submission below propagates the contrast to its
# R scripts via the SLURM environment.
export CONTRAST_NAME MASLD_PROJECT_ROOT="${MASLD_PROJECT_ROOT:-${BASE}}"

# Map labels to on-disk sbatch filenames. Earlier the wrapper invoked names
# like "run_11_perm_fdr.sbatch" that have always carried _{cpu,gpu,io} suffixes
# on disk; submit() silently set JID[label]="" and the pipeline ran 5/10
# pillars. Default to CPU variants; override via env if needed.
declare -A WRAPPERS=(
  [P1]="run_10_M2_refit_random.sbatch"
  [P2]="${WRAPPER_P2:-run_11_perm_cpu.sbatch}"
  [P2_AGG]="run_11b_perm_aggregate.sbatch"
  [P3]="${WRAPPER_P3:-run_12_selection_bias.sbatch}"
  [P3_AGG]="run_12b_aggregate.sbatch"
  [P4]="${WRAPPER_P4:-run_13_boot_cpu.sbatch}"
  [P4_AGG]="run_13b_boot_aggregate.sbatch"
  [P5]="${WRAPPER_P5:-run_14_calibration.sbatch}"
  [P5_AGG]="run_14b_summary.sbatch"
  [P6]="run_15_cochranQ.sbatch"
  [P6_AGG]="run_15b_aggregate.sbatch"
  [P7]="run_16_threshold.sbatch"
  [P8A]="run_17_male_gene_characterization.sbatch"
  [P8B]="run_07b2_loco_magnitude.sbatch"
  [P9_FINAL]="run_18_consensus_v6.sbatch"
)

# Pre-flight: refuse to run if any required sbatch is missing.
missing_wrappers=()
for k in "${!WRAPPERS[@]}"; do
  if [[ ! -e "${SCRIPT_DIR}/${WRAPPERS[$k]}" ]]; then
    missing_wrappers+=("${k}=${WRAPPERS[$k]}")
  fi
done
if [[ ${#missing_wrappers[@]} -gt 0 ]]; then
  echo "FATAL: missing sbatch wrappers (label=filename):"
  printf '  %s\n' "${missing_wrappers[@]}"
  echo "Fix WRAPPERS map or restore the missing files before submitting."
  exit 2
fi

cd "${SCRIPT_DIR}"

# Defaults
FROM_STEP=1
DRY_RUN=0
N_PERM=200
N_BOOT=100
N_SEEDS=50
GAP_FILL=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --from-step) FROM_STEP="$2"; shift 2;;
    --dry-run)   DRY_RUN=1; shift;;
    --n-perm)    N_PERM="$2"; shift 2;;
    --n-boot)    N_BOOT="$2"; shift 2;;
    --n-seeds)   N_SEEDS="$2"; shift 2;;
    --gap-fill)  GAP_FILL=1; shift;;
    -h|--help)   sed -n '1,30p' "$0"; exit 0;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

echo "=== sex_v6 pipeline driver ==="
echo "  from_step=${FROM_STEP}  dry_run=${DRY_RUN}  n_perm=${N_PERM}"
echo "  n_boot=${N_BOOT}  n_seeds=${N_SEEDS}  gap_fill=${GAP_FILL}"
echo "  base=${BASE}"
echo "  sexv3_out=${SEXV3_OUT}"
echo

# ---------------------------------------------------------------------------
# Input file validation
# ---------------------------------------------------------------------------
# Three intermediate RDS files live under sex_v3/intermediates/, and the
# dream-mega-analysis CSV lives at the integration top level. Prior wrapper
# set INPUTS_DIR one level up so all 3 RDS appeared MISSING and the manifest
# emitted "Status: ok / inputs_missing: 3".
# Per-contrast Tier-1 disease-DEG universe (the CSV the pillars read as their
# conditional-filter input). Each contrast has its own disease-side dream output.
case "${CONTRAST_NAME}" in
  disease_vs_ctrl) TIER1_CSV="${INTEGR_DIR}/dream_results_ashr.csv" ;;
  mash_vs_masl)    TIER1_CSV="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/mash_vs_masl_dream.csv" ;;
  masl_vs_ctrl)    TIER1_CSV="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/masl_vs_healthy_dream.csv" ;;
  mash_vs_ctrl)    TIER1_CSV="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/mash_vs_healthy_dream.csv" ;;
  *) echo "FATAL: unknown CONTRAST_NAME='${CONTRAST_NAME}'"; exit 4 ;;
esac
INPUTS=(
  "${INTERMED_DIR}/sex_v3_input.rds"
  "${INTERMED_DIR}/sva_factors.rds"
  "${INTERMED_DIR}/age_mi.rds"
  "${TIER1_CSV}"
)
echo "--- Validating inputs ---"
missing=0
for f in "${INPUTS[@]}"; do
  if [[ -e "${f}" ]]; then
    sz=$(stat -c %s "${f}" 2>/dev/null || echo "?")
    md5=$(md5sum "${f}" 2>/dev/null | awk '{print $1}' || echo "?")
    echo "  OK  ${f}  (${sz} bytes, md5=${md5:0:12})"
  else
    echo "  MISSING ${f}"
    missing=$((missing + 1))
  fi
done
if [[ $missing -gt 0 ]]; then
  echo "FATAL: ${missing} input(s) missing. Refusing to submit; fix inputs first."
  exit 3
fi
echo

mkdir -p "${SCRIPT_DIR}/logs" "${SEXV3_OUT}"

declare -A JID

submit() {
  # submit <label> <sbatch_script> [extra_args...]
  local label="$1"; shift
  local script="$1"; shift
  if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY: would submit ${label}  (${script}) $*"
    JID[$label]="DRY-${label}"
    return
  fi
  if [[ ! -e "${script}" ]]; then
    echo "MISSING sbatch ${label}: ${script}  (skipping)"
    JID[$label]=""
    return
  fi
  local out
  # Always export CONTRAST_NAME so per-pillar R scripts pick up the contrast
  # via contrast_paths(). --export=ALL preserves the current shell environment
  # plus the explicit variables.
  out=$(sbatch --parsable --export=ALL,CONTRAST_NAME="${CONTRAST_NAME}" "$@" "${script}")
  JID[$label]="${out}"
  echo "submitted ${label} jid=${out}"
}

# Build pipeline plan
SUBMIT_P1=0; SUBMIT_P3=0; SUBMIT_P5=0; SUBMIT_P6=0
SUBMIT_P2=0; SUBMIT_P4=0; SUBMIT_P7=0
SUBMIT_8A=0; SUBMIT_8B=0; SUBMIT_FINAL=0

if [[ $FROM_STEP -le 1 ]]; then SUBMIT_P1=1; fi
if [[ $FROM_STEP -le 2 ]]; then SUBMIT_P3=1; fi
if [[ $FROM_STEP -le 3 ]]; then SUBMIT_P5=1; fi
if [[ $FROM_STEP -le 4 ]]; then SUBMIT_P6=1; fi
if [[ $FROM_STEP -le 5 ]]; then SUBMIT_P2=1; fi
if [[ $FROM_STEP -le 6 ]]; then SUBMIT_P4=1; fi
if [[ $FROM_STEP -le 7 ]]; then SUBMIT_P7=1; fi
if [[ $FROM_STEP -le 8 ]]; then SUBMIT_8A=1; SUBMIT_8B=1; fi
if [[ $FROM_STEP -le 9 ]]; then SUBMIT_FINAL=1; fi

echo "--- Submitting jobs ---"

# Helper: build "afterok:JID1:JID2:..." from one or more labels (skip empty / DRY).
dep_afterok() {
  local kind="$1"; shift
  local jids=""
  for k in "$@"; do
    local v="${JID[$k]:-}"
    if [[ -n "$v" && "${v:0:4}" != "DRY-" ]]; then
      jids+="${v}:"
    fi
  done
  jids="${jids%:}"
  if [[ -n "$jids" ]]; then
    echo "--dependency=${kind}:${jids}"
  fi
}

[[ $SUBMIT_P1 -eq 1 ]] && submit P1 "${WRAPPERS[P1]}"
[[ $SUBMIT_P3 -eq 1 ]] && submit P3 "${WRAPPERS[P3]}"
[[ $SUBMIT_P5 -eq 1 ]] && submit P5 "${WRAPPERS[P5]}" --export=ALL,N_SEEDS=${N_SEEDS}
[[ $SUBMIT_P6 -eq 1 ]] && submit P6 "${WRAPPERS[P6]}"

P1_DEP=$(dep_afterok afterok P1)
[[ -z "$P1_DEP" && $DRY_RUN -eq 1 ]] && P1_DEP="--dependency=afterok:DRY-P1"

[[ $SUBMIT_P2 -eq 1 ]] && submit P2 "${WRAPPERS[P2]}" ${P1_DEP} --export=ALL,N_PERM=${N_PERM}
[[ $SUBMIT_P4 -eq 1 ]] && submit P4 "${WRAPPERS[P4]}" ${P1_DEP} --export=ALL,N_BOOT=${N_BOOT}
[[ $SUBMIT_P7 -eq 1 ]] && submit P7 "${WRAPPERS[P7]}" ${P1_DEP}

# Aggregators (B0-7) — wire afterok deps so each aggregator runs only after
# its worker array completes successfully. v6 results were previously rebuilt
# by hand-submitted JIDs because these were not in the dep graph.
P2_DEP=$(dep_afterok afterok P2)
P3_DEP=$(dep_afterok afterok P3)
P4_DEP=$(dep_afterok afterok P4)
P5_DEP=$(dep_afterok afterok P5)
P6_DEP=$(dep_afterok afterok P6)

[[ $SUBMIT_P2 -eq 1 ]] && submit P2_AGG "${WRAPPERS[P2_AGG]}" ${P2_DEP}
[[ $SUBMIT_P3 -eq 1 ]] && submit P3_AGG "${WRAPPERS[P3_AGG]}" ${P3_DEP}
[[ $SUBMIT_P4 -eq 1 ]] && submit P4_AGG "${WRAPPERS[P4_AGG]}" ${P4_DEP}
[[ $SUBMIT_P5 -eq 1 ]] && submit P5_AGG "${WRAPPERS[P5_AGG]}" ${P5_DEP}
[[ $SUBMIT_P6 -eq 1 ]] && submit P6_AGG "${WRAPPERS[P6_AGG]}" ${P6_DEP}

P6_AGG_DEP=$(dep_afterok afterany P6_AGG)
[[ -z "$P6_AGG_DEP" && $DRY_RUN -eq 1 ]] && P6_AGG_DEP="--dependency=afterany:DRY-P6_AGG"

[[ $SUBMIT_8A -eq 1 ]] && submit P8A "${WRAPPERS[P8A]}" ${P6_AGG_DEP}
[[ $SUBMIT_8B -eq 1 ]] && submit P8B "${WRAPPERS[P8B]}"

# Final consensus: depends on all aggregators (canonical inputs into 18).
FINAL_DEP=$(dep_afterok afterany P1 P2_AGG P3_AGG P4_AGG P5_AGG P6_AGG P7 P8A P8B)
[[ $SUBMIT_FINAL -eq 1 ]] && submit P9_FINAL "${WRAPPERS[P9_FINAL]}" ${FINAL_DEP}

# ---------------------------------------------------------------------------
# Initial manifest scaffold
# ---------------------------------------------------------------------------
{
  echo "{"
  echo "  \"created_utc\": \"$(date -u +%Y-%m-%dT%H:%M:%SZ)\","
  echo "  \"from_step\": ${FROM_STEP},"
  echo "  \"dry_run\": ${DRY_RUN},"
  echo "  \"n_perm\": ${N_PERM},"
  echo "  \"n_boot\": ${N_BOOT},"
  echo "  \"n_seeds\": ${N_SEEDS},"
  echo "  \"inputs_missing\": ${missing},"
  echo "  \"job_ids\": {"
  first=1
  for k in P1 P3 P5 P6 P2 P4 P7 P8A P8B P2_AGG P3_AGG P4_AGG P5_AGG P6_AGG P9_FINAL; do
    v="${JID[$k]:-null}"
    if [[ "$v" == "null" || -z "$v" ]]; then v="null"; else v="\"$v\""; fi
    if [[ $first -eq 1 ]]; then first=0; else echo ","; fi
    printf "    \"%s\": %s" "$k" "$v"
  done
  echo
  echo "  },"
  echo "  \"status\": \"submitted\""
  echo "}"
} > "${MANIFEST}"

echo
echo "--- Submission summary ---"
for k in P1 P3 P5 P6 P2 P4 P7 P8A P8B P2_AGG P3_AGG P4_AGG P5_AGG P6_AGG P9_FINAL; do
  echo "  ${k}: ${JID[$k]:-(not submitted)}"
done
echo
echo "Manifest scaffold: ${MANIFEST}"
echo "Pipeline driver done."
