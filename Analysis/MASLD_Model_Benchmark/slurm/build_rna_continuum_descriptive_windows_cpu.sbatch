#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
HAC="${PROJECT}/RNA-seq/results/histology_anchored_continuum"

# Job 1 output. The clinical values this layer reports are the ones re-derived and
# verified there, not the manifest's unverified copies.
PROVENANCE="${ROOT}/executions/rna-continuum-clinical-provenance-21113968"
PROVENANCE_SHA256=4d5e2349f85602d03ec0c82ab76a0810296a8cae3480b263a322b0e2000931aa

CANDIDATE="${HAC}/candidates/hac-continuum-20260818T024923Z"
LAYERS="${HAC}/molecular_layers/hac-molecular-layers-20260818T173348Z"
SYSTEMS="${HAC}/molecular_systems_inference/hac-system-inference-20260824T232135Z"

AXIS_SCORES="${CANDIDATE}/projection/participant_scores.tsv"
PROGRAM_TESTABILITY="${CANDIDATE}/programs/program_testability.tsv"
HOTSPOT_SCORES="${LAYERS}/programs/hotspot/hotspot_participant_scores.tsv.gz"
FROZEN_WINDOWS="${LAYERS}/programs/hotspot/hotspot_fixed_windows.tsv.gz"
SYSTEM_SCORES="${SYSTEMS}/tables/participant_system_scores.tsv.gz"

TARGET="${ROOT}/executions/rna-continuum-descriptive-windows-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite descriptive window artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
mkdir -p "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"

cat > "${STAGE}/validation/bound_fixtures.sha256" <<SUMS
${PROVENANCE_SHA256}  ${PROVENANCE}/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

# Record the upstream score output files as observed. They live outside the benchmark and
# are frozen candidate trees, so their digests are captured rather than pre-pinned.
sha256sum "${AXIS_SCORES}" "${PROGRAM_TESTABILITY}" "${HOTSPOT_SCORES}" \
          "${FROZEN_WINDOWS}" "${SYSTEM_SCORES}" \
  > "${STAGE}/validation/upstream_score_inputs.sha256"

python3 -m unittest -v tests.unit.test_build_rna_continuum_descriptive_windows \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

python3 scripts/build_rna_continuum_descriptive_windows.py \
  --provenance-artifact "${PROVENANCE}" \
  --provenance-sha256 "${PROVENANCE_SHA256}" \
  --axis-scores "${AXIS_SCORES}" \
  --hotspot-scores "${HOTSPOT_SCORES}" \
  --system-scores "${SYSTEM_SCORES}" \
  --program-testability "${PROGRAM_TESTABILITY}" \
  --frozen-hotspot-windows "${FROZEN_WINDOWS}" \
  --output "${STAGE}" \
  > "${STAGE}/validation/windows_stdout.json"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum scripts/build_rna_continuum_descriptive_windows.py \
          tests/unit/test_build_rna_continuum_descriptive_windows.py \
          slurm/build_rna_continuum_descriptive_windows_cpu.sbatch > "${STAGE}/source.sha256"

python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "provenance" / "audit.json").read_text())

if audit["donor_level_claim_made"] is not False:
    raise SystemExit("this artifact may not make a donor-level claim")
if audit["predictive_arm_built"] is not False:
    raise SystemExit("a predictive continuum arm is barred by the frozen contract")
if audit["inference_use"] is not False:
    raise SystemExit("sliding windows are a reporting view and carry inference_use false")
if audit["arms_pooled"] is not False:
    raise SystemExit("the quarantined arm may never be pooled with the primary arm")
if audit["reproduction_control"]["n_mismatched"] != 0:
    raise SystemExit("recomputed windows disagree with the frozen upstream table")
if audit["reproduction_control"]["n_checked"] == 0:
    raise SystemExit("the reproduction control checked nothing")
if audit["features_expected"] != audit["features_emitted"]:
    raise SystemExit(
        f"feature shortfall: expected {audit['features_expected']} features, "
        f"emitted {audit['features_emitted']}"
    )
for check in audit["shape_checks"]:
    if check["rows_expected"] != check["rows_emitted"]:
        raise SystemExit(
            f"{check['arm']}: expected {check['rows_expected']} window rows, "
            f"emitted {check['rows_emitted']}"
        )

freeze_tree(root, {
    "artifact_class": "rna_continuum_descriptive_window_layer",
    "workstream": "H_continuum_rna_descriptive",
    "unit_of_observation": "sequencing_sample",
    "donor_level_claim_made": False,
    "predictive_arm_built": False,
    "inference_use": False,
    "arms_pooled": False,
    "axis_id": audit["axis_id"],
    "n_window_rows_primary": audit["n_window_rows_primary"],
    "n_window_rows_quarantined": audit["n_window_rows_quarantined"],
    "reproduction_control_checked": audit["reproduction_control"]["n_checked"],
    "provenance_artifact_sha256": audit["provenance_artifact_sha256"],
    "features_expected": audit["features_expected"],
    "features_emitted": audit["features_emitted"],
    "supersedes": audit["supersedes"],
    "supersedes_reason": audit["supersedes_reason"],
    "release_state": "candidate_only",
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Descriptive window layer: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
