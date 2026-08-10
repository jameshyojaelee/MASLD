#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 RUN_ROOT" >&2
  exit 64
fi
RUN_ROOT=$(realpath "$1")
[[ -f "$RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 64; }
[[ -f "$RUN_ROOT/RECOUNT_COMPLETE" ]] || { echo "Recount is not complete" >&2; exit 64; }
[[ -f "$RUN_ROOT/BASELINE_FROZEN.json" ]] || { echo "Baseline is not frozen" >&2; exit 64; }
[[ ! -e "$RUN_ROOT/RECOUNT_IMPORT_IN_PROGRESS" && ! -e "$RUN_ROOT/RECOUNT_IMPORT_FAILED.json" ]] || {
  echo "Imported recount is incomplete or failed" >&2
  exit 64
}
[[ ! -e "$RUN_ROOT/ANALYSIS_COMPLETE.json" ]] || { echo "Refusing to reuse a completed analysis candidate" >&2; exit 64; }
[[ ! -e "$RUN_ROOT/ANALYSIS_INVALID.json" ]] || { echo "Refusing to reuse an INVALID analysis candidate" >&2; exit 64; }

SNAP="$RUN_ROOT/source_snapshot"
PROJECT_ROOT="$SNAP"
# Harness code defaults to the frozen snapshot. BG001_SCRIPT_DIR allows an
# explicitly authorized deviation (e.g. re-running the analysis after a guard
# fix that the already-frozen snapshot predates) without touching the snapshot,
# so verify_run_contract.py still verifies every frozen byte. ANALYSIS INPUTS
# are unaffected: INTEGRATION, CONFIG, METADATA and the frozen sets below all
# remain pinned to "$SNAP" regardless of this override.
SCRIPT_DIR="${BG001_SCRIPT_DIR:-$SNAP/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation}"
# Producers default to the frozen snapshot. BG001_INTEGRATION_DIR permits an
# explicitly authorized deviation for a producer fix the snapshot predates. Use
# only with a recorded per-file diff; DATA inputs (CONFIG, METADATA, frozen_sets)
# stay pinned to "$SNAP" regardless.
INTEGRATION="${BG001_INTEGRATION_DIR:-$SNAP/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts}"
CONFIG="$SNAP/config/human_datasets.yaml"
METADATA="$SNAP/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
OVERRIDES="$RUN_ROOT/manifests/count_overrides.tsv"
R0_OVERRIDES="$RUN_ROOT/manifests/read_count_overrides.tsv"
LOCKED_QC="$RUN_ROOT/frozen_sets/locked_sample_qc_report.csv"
LOCKED_META="$RUN_ROOT/frozen_sets/locked_meta_matched.rds"
LOCKED_DGE="$RUN_ROOT/frozen_sets/locked_merged_dge.rds"
REFERENCE_DEG="$RUN_ROOT/frozen_sets/canonical_deg_results.csv"
GENE_METADATA="$SNAP/data/gencode_v49_gene_metadata.tsv.gz"

for required in "$OVERRIDES" "$R0_OVERRIDES" "$CONFIG" "$METADATA" "$LOCKED_QC" "$LOCKED_META" "$LOCKED_DGE" "$REFERENCE_DEG" "$GENE_METADATA"; do
  [[ -f "$required" ]] || { echo "Missing required input: $required" >&2; exit 64; }
done
for arm in R0 F_locked F_legacy F_five; do
  arm_root="$RUN_ROOT/arms/$arm"
  mapfile -d '' -t existing_arm_files < <(find "$arm_root" -type f ! -name '.bg001_candidate_root' -print0)
  [[ ${#existing_arm_files[@]} -eq 0 ]] || {
    echo "Refusing nonempty or previously attempted arm $arm: ${existing_arm_files[0]}" >&2
    exit 64
  }
done

eval "$(micromamba shell hook --shell=bash)"
set +u
micromamba activate rnaseq
set -u
python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$RUN_ROOT" --require-bam-frozen
python3 "$SCRIPT_DIR/verify_analysis_contract.py" \
  --run-root "$RUN_ROOT" --require-baseline-frozen
python3 "$SCRIPT_DIR/analysis_runtime_contract.py" verify --run-root "$RUN_ROOT"
R_SCRIPT=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["Rscript"]["path"])' \
  "$RUN_ROOT/contract/analysis_runtime_contract.json")
SEAL_BOUND_ARGS=(
  --bound-input "$RUN_ROOT/contract/run_contract.json"
  --bound-input "$RUN_ROOT/contract/source_manifest.tsv"
  --bound-input "$RUN_ROOT/contract/analysis_runtime_contract.json"
  --bound-input "$RUN_ROOT/contract/analysis_environment.explicit.txt"
  --bound-input "$RUN_ROOT/contract/analysis_environment.conda.json"
  --bound-input "$RUN_ROOT/contract/analysis_r_dependency_files.tsv"
  --bound-input "$RUN_ROOT/contract/analysis_native_dependencies.tsv"
  --bound-input "$RUN_ROOT/BASELINE_FROZEN.json"
)

verify_runtime() {
  python3 "$SCRIPT_DIR/analysis_runtime_contract.py" verify --run-root "$RUN_ROOT"
}

run_pinned_r() {
  verify_runtime
  set +e
  "$R_SCRIPT" "$@"
  local status=$?
  set -e
  verify_runtime
  return "$status"
}

run_pinned_r_env() {
  local -a environment=()
  while [[ $# -gt 0 && $1 != -- ]]; do
    environment+=("$1")
    shift
  done
  [[ $# -gt 0 && $1 == -- ]] || { echo "run_pinned_r_env requires --" >&2; return 64; }
  shift
  verify_runtime
  set +e
  env "${environment[@]}" "$R_SCRIPT" "$@"
  local status=$?
  set -e
  verify_runtime
  return "$status"
}

verify_recount() {
  python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$RUN_ROOT" --require-bam-frozen
  python3 "$SCRIPT_DIR/verify_analysis_contract.py" \
    --run-root "$RUN_ROOT" --require-baseline-frozen
  verify_runtime
  if [[ -e "$RUN_ROOT/RECOUNT_IMPORTED" ]]; then
    python3 "$SCRIPT_DIR/import_recount.py" --verify-import --destination-run "$RUN_ROOT"
  else
    python3 "$SCRIPT_DIR/verify_recount_artifacts.py" --run-root "$RUN_ROOT"
  fi
}

run_qc() {
  local arm=$1 mode=$2 overrides=${3:-}
  local arm_root="$RUN_ROOT/arms/$arm"
  echo "[$arm] 01 QC mode=$mode"
  run_pinned_r_env \
    MASLD_PROJECT_ROOT="$PROJECT_ROOT" \
    MASLD_RUN_ROOT="$arm_root" \
    MASLD_CONFIG_PATH="$CONFIG" \
    MASLD_METADATA_PATH="$METADATA" \
    MASLD_QC_MODE="$mode" \
    MASLD_COUNT_OVERRIDES="$overrides" \
    MASLD_LOCKED_QC="$LOCKED_QC" \
    MASLD_LOCKED_META="$LOCKED_META" \
    -- --vanilla "$INTEGRATION/01_sample_qc.R" \
    >"$arm_root/logs/01_sample_qc.log" 2>&1
}

run_dge() {
  local arm=$1 input_arm=$2 scope=$3 gene_mode=$4 counts_mode=$5 qc_mode=$6
  local arm_root="$RUN_ROOT/arms/$arm"
  local input_root="$RUN_ROOT/arms/$input_arm"
  echo "[$arm] 03 filter_scope=$scope gene_mode=$gene_mode input=$input_arm"
  run_pinned_r_env \
    MASLD_PROJECT_ROOT="$PROJECT_ROOT" \
    MASLD_RUN_ROOT="$arm_root" \
    MASLD_INPUT_RUN_ROOT="$input_root" \
    MASLD_CONFIG_PATH="$CONFIG" \
    MASLD_FILTER_SCOPE="$scope" \
    MASLD_GENE_MODE="$gene_mode" \
    MASLD_COUNTS_MODE="$counts_mode" \
    MASLD_QC_MODE="$qc_mode" \
    MASLD_LOCKED_DGE="$LOCKED_DGE" \
    -- --vanilla "$INTEGRATION/03_integrate_counts.R" \
    >"$arm_root/logs/03_integrate_counts.log" 2>&1
}

run_fit() {
  local arm=$1 meta_arm=${2:-$1}
  local arm_root="$RUN_ROOT/arms/$arm"
  local meta_root="$RUN_ROOT/arms/$meta_arm"
  echo "[$arm] 05h C2 LVQW/TREAT"
  run_pinned_r_env \
    MASLD_PROJECT_ROOT="$PROJECT_ROOT" \
    MASLD_RUN_ROOT="$arm_root" \
    MASLD_CONFIG_PATH="$CONFIG" \
    MASLD_DGE_INPUT="$arm_root/results/integration/merged_dge.rds" \
    MASLD_META_INPUT="$meta_root/results/integration/meta_matched.rds" \
    MASLD_DEG_OUTPUT="$arm_root/results/integration/deg_results.csv" \
    MASLD_SANITY_REFERENCE="$REFERENCE_DEG" \
    MASLD_GENE_METADATA_PATH="$GENE_METADATA" \
    CANONICAL_TREAT_LFC=0.25 \
    -- --vanilla "$INTEGRATION/05h_limma_voom_qw_canonical.R" \
    >"$arm_root/logs/05h_limma_voom_qw.log" 2>&1
}

freeze_arm() {
  local arm=$1
  local arm_root="$RUN_ROOT/arms/$arm"
  local manifest="$arm_root/provenance/artifact_manifest.tsv"
  local marker="$arm_root/ARM_COMPLETE.json"
  [[ ! -e "$manifest" && ! -e "$marker" ]] || {
    echo "Refusing to refreeze completed or partially frozen arm $arm" >&2
    exit 64
  }
  verify_runtime
  python3 "$SCRIPT_DIR/seal_analysis_artifacts.py" seal \
    --run-root "$RUN_ROOT" \
    --artifact-root "$arm_root" \
    --manifest "$manifest" \
    --marker "$marker" \
    --label "$arm" \
    "${SEAL_BOUND_ARGS[@]}"
}

verify_arm() {
  local arm=$1
  local arm_root="$RUN_ROOT/arms/$arm"
  python3 "$SCRIPT_DIR/seal_analysis_artifacts.py" verify \
    --run-root "$RUN_ROOT" \
    --artifact-root "$arm_root" \
    --manifest "$arm_root/provenance/artifact_manifest.tsv" \
    --marker "$arm_root/ARM_COMPLETE.json" \
    --label "$arm"
}

verify_recount
run_qc R0 recompute "$R0_OVERRIDES"
run_dge R0 R0 legacy_all native read recompute
run_fit R0
run_pinned_r --vanilla "$SCRIPT_DIR/validate_arm.R" "$RUN_ROOT" R0
set +e
python3 "$SCRIPT_DIR/check_r0_reproduction.py" --run-root "$RUN_ROOT"
r0_status=$?
set -e
if [[ $r0_status -ne 0 ]]; then
  python3 "$SCRIPT_DIR/seal_analysis_invalid.py" seal \
    --run-root "$RUN_ROOT" --kind r0
  python3 "$SCRIPT_DIR/seal_analysis_invalid.py" verify \
    --run-root "$RUN_ROOT" --kind r0
  exit "$r0_status"
fi
run_pinned_r --vanilla "$SCRIPT_DIR/check_locked_substrate.R" "$RUN_ROOT" R0
freeze_arm R0

verify_recount
run_qc F_locked locked "$OVERRIDES"
run_pinned_r --vanilla "$SCRIPT_DIR/check_locked_substrate.R" "$RUN_ROOT" F_locked
run_dge F_locked F_locked legacy_all locked fragment locked
run_fit F_locked
run_pinned_r --vanilla "$SCRIPT_DIR/validate_arm.R" "$RUN_ROOT" F_locked
freeze_arm F_locked

verify_recount
run_qc F_legacy recompute "$OVERRIDES"
run_dge F_legacy F_legacy legacy_all native fragment recompute
run_fit F_legacy
run_pinned_r --vanilla "$SCRIPT_DIR/validate_arm.R" "$RUN_ROOT" F_legacy
freeze_arm F_legacy

# F_five reuses the accepted F_legacy QC/meta/raw-count objects, but copies them
# into its own arm so every candidate is self-contained and hashable.
verify_recount
cp -a "$RUN_ROOT/arms/F_legacy/qc/." "$RUN_ROOT/arms/F_five/qc/"
cp -p "$RUN_ROOT/arms/F_legacy/results/integration/merged_counts_raw.rds" "$RUN_ROOT/arms/F_five/results/integration/merged_counts_raw.rds"
cp -p "$RUN_ROOT/arms/F_legacy/results/integration/meta_matched.rds" "$RUN_ROOT/arms/F_five/results/integration/meta_matched.rds"
cp -p "$RUN_ROOT/arms/F_legacy/provenance/effective_count_sources.tsv" "$RUN_ROOT/arms/F_five/provenance/effective_count_sources.tsv"
run_dge F_five F_five canonical_five native fragment reuse_F_legacy
run_fit F_five
run_pinned_r --vanilla "$SCRIPT_DIR/validate_arm.R" "$RUN_ROOT" F_five
freeze_arm F_five

for arm in R0 F_locked F_legacy F_five; do verify_arm "$arm"; done

verify_recount
run_pinned_r --vanilla "$SCRIPT_DIR/compare_scientific_gates.R" "$SNAP" "$RUN_ROOT"
python3 "$SCRIPT_DIR/generate_dependency_manifest.py" --project-root "$SNAP" --run-root "$RUN_ROOT"
verify_runtime
COMPARISON_MANIFEST="$RUN_ROOT/comparisons/artifact_manifest.tsv"
ANALYSIS_MARKER="$RUN_ROOT/ANALYSIS_COMPLETE.json"
FINAL_BOUND_ARGS=("${SEAL_BOUND_ARGS[@]}" --bound-input "$RUN_ROOT/RECOUNT_COMPLETE")
for arm in R0 F_locked F_legacy F_five; do
  FINAL_BOUND_ARGS+=(
    --bound-input "$RUN_ROOT/arms/$arm/provenance/artifact_manifest.tsv"
    --bound-input "$RUN_ROOT/arms/$arm/ARM_COMPLETE.json"
  )
done
[[ ! -e "$COMPARISON_MANIFEST" && ! -e "$ANALYSIS_MARKER" ]] || {
  echo "Refusing an existing comparison artifact manifest" >&2
  exit 64
}
python3 "$SCRIPT_DIR/seal_analysis_artifacts.py" seal \
  --run-root "$RUN_ROOT" \
  --artifact-root "$RUN_ROOT/comparisons" \
  --manifest "$COMPARISON_MANIFEST" \
  --marker "$ANALYSIS_MARKER" \
  --label comparisons \
  "${FINAL_BOUND_ARGS[@]}"
for arm in R0 F_locked F_legacy F_five; do verify_arm "$arm"; done
verify_recount
python3 "$SCRIPT_DIR/seal_analysis_artifacts.py" verify \
  --run-root "$RUN_ROOT" \
  --artifact-root "$RUN_ROOT/comparisons" \
  --manifest "$COMPARISON_MANIFEST" \
  --marker "$ANALYSIS_MARKER" \
  --label comparisons
verify_runtime
