#!/usr/bin/env bash
#SBATCH --job-name=model-eval-894
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
DECONV=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution
TASKSPEC=${ROOT}/config/evaluation/bulk_lineage_composition_fibrosis_task.toml
TASKSPEC_SHA256=42f7b71a7f145d10cf2fcddc8b41da015db2f699361a5ad9ce49edc5622ad5b0
GATE=${ROOT}/config/evaluation/bulk_lineage_composition_fibrosis_promotion_gate.json
GATE_SHA256=e9cb9f9d11d1023450b5a50115f0b1b61b326caf99a84c12e1abee88b1cc05d9
ACTIVATION=${ROOT}/config/evaluation/bulk_existing_source_pool_composition_activation.json
ACTIVATION_SHA256=6fb3cccd8c79d3e445725329086f8adfe1904753d2fd3c5c975be0ec06037de1
TARGET=${ROOT}/executions/model-eval-894-${SLURM_JOB_ID}-lineage-composition
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the lineage-composition evaluation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${TASKSPEC_SHA256}" "${TASKSPEC}"     | sha256sum --check --strict
printf '%s  %s\n' "${GATE_SHA256}" "${GATE}"             | sha256sum --check --strict
printf '%s  %s\n' "${ACTIVATION_SHA256}" "${ACTIVATION}" | sha256sum --check --strict
bash -n slurm/run_bulk_lineage_composition_fibrosis_cpu.sbatch

# The read-only deconvolution tree must not be written.
find "${DECONV}/bulk" "${DECONV}/results" -maxdepth 2 -type f \
  \( -name '*_metadata.tsv' -o -name '*_bayesprism_proportions.tsv' \) \
  -printf '%p %s\n' | sort > "${STAGE}/source_tree_before.txt"

# MuSiC must never be an input.
if grep -rn 'music' scripts/analyze_bulk_lineage_composition_fibrosis.py \
     scripts/freeze_lineage_composition_substrate.py | grep -vi 'music_used\|never MuSiC\|MuSiC output'; then
  printf 'A MuSiC path unexpectedly appears in the lane scripts.\n' >&2
  exit 1
fi

"${PYTHON}" -m unittest -v tests.unit.test_bulk_lineage_composition_fibrosis \
  > "${STAGE}/unit_tests.txt" 2>&1

"${PYTHON}" scripts/freeze_lineage_composition_substrate.py \
  --deconvolution-root "${DECONV}" \
  --output "${STAGE}/substrate" \
  > "${STAGE}/substrate_stdout.json" 2> "${STAGE}/substrate_stderr.txt"

"${PYTHON}" scripts/analyze_bulk_lineage_composition_fibrosis.py \
  --fixture "${STAGE}/substrate/fixture" \
  --gate "${GATE}" \
  --gate-sha256 "${GATE_SHA256}" \
  --output "${STAGE}/evaluation" \
  --seed 20260826 \
  --bootstrap-replicates 10000 \
  --permutation-replicates 10000 \
  > "${STAGE}/evaluation_stdout.json" 2> "${STAGE}/evaluation_stderr.txt"

find "${DECONV}/bulk" "${DECONV}/results" -maxdepth 2 -type f \
  \( -name '*_metadata.tsv' -o -name '*_bayesprism_proportions.tsv' \) \
  -printf '%p %s\n' | sort > "${STAGE}/source_tree_after.txt"
if ! cmp -s "${STAGE}/source_tree_before.txt" "${STAGE}/source_tree_after.txt"; then
  printf 'The read-only deconvolution tree changed during the run.\n' >&2
  exit 1
fi

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/freeze_lineage_composition_substrate.py \
  scripts/analyze_bulk_lineage_composition_fibrosis.py \
  tests/unit/test_bulk_lineage_composition_fibrosis.py \
  slurm/run_bulk_lineage_composition_fibrosis_cpu.sbatch \
  "${TASKSPEC}" "${GATE}" "${ACTIVATION}" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
substrate = root / "substrate"
evaluation = root / "evaluation"
sub = json.loads((substrate / "substrate_receipt.json").read_text(encoding="utf-8"))
rec = json.loads((evaluation / "receipt.json").read_text(encoding="utf-8"))
verdict = json.loads((evaluation / "promotion_gate_verdict.json").read_text(encoding="utf-8"))

if (
    sub.get("status") != "pass_substrate_frozen_hashed_copy"
    or sub.get("lineages") != 16
    or sub.get("primary_total_samples") != 350
    or sub.get("control_arms_excluded_total") != 39
    or sub.get("donor_key_exists") is not False
    or sub.get("music_used") is not False
    or rec.get("status") != "passed_external_development_scoring"
    or rec.get("lineages") != 16
    or rec.get("detection_floor") != 1e-4
    or rec.get("donor_count_claimed") is not False
    or rec.get("cohort_level_bootstrap_performed") is not False
    or rec.get("between_cohort_heterogeneity_estimated") is not False
    or rec.get("combination") != "macro_average_across_cohort_families"
    or rec.get("champion_claim_eligible") is not False
    or rec["two_cohort_co_primary"]["is_not_a_fresh_frozen_estimate"] is not True
    or verdict.get("verdict") not in {"PASS", "FAIL"}
    or verdict.get("conditions_applicable") < 1
):
    raise SystemExit("lineage-composition verification differs")
for directory, cls in (
    (substrate, "lineage_composition_substrate"),
    (evaluation, "lineage_composition_evaluation"),
):
    freeze_tree(directory, {
        "artifact_class": cls,
        "lineages": 16,
        "primary_total_samples": 350,
        "promotion_gate_verdict": verdict["verdict"],
        "external_development_only": True,
        "champion_eligible": False,
        "status": "passed",
    })
    verify_frozen_tree(directory)
freeze_tree(root, {
    "artifact_class": "bulk_lineage_composition_fibrosis_campaign",
    "promotion_gate_verdict": verdict["verdict"],
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps(verdict, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/substrate/ARTIFACTS.json" "${TARGET}/evaluation/ARTIFACTS.json"
