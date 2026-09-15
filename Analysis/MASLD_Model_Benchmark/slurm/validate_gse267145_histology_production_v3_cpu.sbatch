#!/usr/bin/env bash
#SBATCH --job-name=model-check-072
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export OPENBLAS_NUM_THREADS=3

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds
FIT_VIEWS=${ROOT}/executions/model-training-067-21080776/fit_views
SURFACE=${ROOT}/config/evaluation/gse267145_histology_benchmark_surface.json
CONTRACT=${ROOT}/config/evaluation/gse267145_histology_production_v3.json
FITTER=${ROOT}/scripts/fit_gse267145_histology_baselines_v3.py
V2=${ROOT}/executions/model-training-069-production-v2
TARGET=${ROOT}/executions/model-check-072-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 v3 validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/units"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

for pair in \
  "8be11a742def0f80465bca92fd02baf626d2d207ee30a87cdb25bcf59e47e29f config/evaluation/gse267145_histology_production_v3.json" \
  "0cf28366d86040ce18fec73e3a5cde862ab33dad71cf0719bdd1d15a14508d3e scripts/fit_gse267145_histology_baselines_v3.py" \
  "c75dc8dee28231e1377bb1cf8a0c5d269a29fe0aaa94be71ecaf37c6cfea0240 tests/unit/test_fit_gse267145_histology_baselines_v3.py" \
  "b9832e9ddb1fe794103167db23467099288e3d9c351971716e043d6405d9dd27 scripts/fit_gse267145_histology_baselines.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

"${PYTHON}" -m py_compile \
  scripts/fit_gse267145_histology_baselines.py \
  scripts/fit_gse267145_histology_baselines_v3.py
"${PYTHON}" -m unittest -v \
  tests.unit.test_fit_gse267145_histology_baselines \
  tests.unit.test_fit_gse267145_histology_baselines_v3 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

declare -a units=("1 1701" "1 1721")
declare -a pids=()
for unit in "${units[@]}"; do
  read -r outer seed <<< "${unit}"
  unit_root=${STAGE}/units/outer_${outer}--seed_${seed}
  install -d -m 0750 "${unit_root}"
  /usr/bin/time -f \
    '{"elapsed_seconds": %e, "maximum_resident_kib": %M, "user_seconds": %U, "system_seconds": %S}' \
    -o "${unit_root}/resource.json" \
    "${PYTHON}" "${FITTER}" \
      --molecular "${MOLECULAR}" \
      --molecular-artifacts-sha256 2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456 \
      --folds "${FOLDS}" \
      --folds-artifacts-sha256 9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3 \
      --fit-views "${FIT_VIEWS}" \
      --fit-views-artifacts-sha256 8254274f1b11089656f4b6d5ec446a6ba491e2f28ae8abcc088fc6e78ac75d6a \
      --surface "${SURFACE}" \
      --surface-sha256 23367a28bd16266c79adde854534b124cfad7d2062788625839e78f3267a8c0b \
      --preflight-outer-fold "${outer}" \
      --preflight-seed "${seed}" \
      --output "${unit_root}/fit" \
      > "${unit_root}/stdout.json" 2> "${unit_root}/stderr.txt" &
  pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then status=1; fi
done
if [[ "${status}" -ne 0 ]]; then
  printf 'At least one GSE267145 v3 validation unit failed.\n' >&2
  exit 1
fi

for prediction in "${V2}/units/outer_1/seed_1701/fit/predictions/"*.tsv; do
  cmp --silent \
    "${prediction}" \
    "${STAGE}/units/outer_1--seed_1701/fit/predictions/${prediction##*/}"
done

"${PYTHON}" - "${STAGE}" "${CONTRACT}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
contract = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
expected = {(1, 1701): 0, (1, 1721): 1}
observed = {}
for path in sorted((root / "units").glob("outer_*--seed_*")):
    outer = int(path.name.split("--")[0].removeprefix("outer_"))
    seed = int(path.name.split("--")[1].removeprefix("seed_"))
    receipt = json.loads((path / "fit/receipt.json").read_text(encoding="utf-8"))
    fallback = receipt.get("fibrosis_group3_fallback", {})
    activation_count = fallback.get("activation_count")
    if (
        receipt.get("schema_version")
        != "masld-bench-gse267145-histology-baseline-fit-preflight-v3"
        or receipt.get("revision_id") != contract.get("revision_id")
        or receipt.get("status") != "passed_timing_predictions_unscored"
        or receipt.get("outer_fold") != outer
        or receipt.get("seed") != seed
        or receipt.get("exact_production_grids_used") is not True
        or receipt.get("prediction_files") != 11
        or receipt.get("primary_endpoint_code_path_changed") is not False
        or receipt.get("other_secondary_endpoint_code_paths_changed") is not False
        or receipt.get("outer_test_outcomes_read") is not False
        or receipt.get("outer_test_metrics_calculated") is not False
        or fallback.get("fallback_id")
        != "training_only_fibrosis_group3_class_prior_v1"
        or fallback.get("eligible_endpoint") != "fibrosis_group3"
        or fallback.get("eligible_endpoint_role") != "secondary_only"
        or fallback.get("outer_test_outcomes_used") is not False
        or fallback.get("metrics_calculated") is not False
    ):
        raise SystemExit(f"v3 validation receipt differs: {path}")
    minimum = expected[(outer, seed)]
    if seed == 1701 and activation_count != 0:
        raise SystemExit("v3 changed the previously successful 1/1701 unit")
    if seed == 1721 and (not isinstance(activation_count, int) or activation_count < minimum):
        raise SystemExit("v3 fallback did not activate for 1/1721")
    selections = json.loads(
        (
            path
            / f"fit/fit_receipts/outer_{outer}/seed_{seed}/selection_receipts.json"
        ).read_text(encoding="utf-8")
    )
    activated = sorted(fallback.get("activated_model_ids", []))
    for model_id in activated:
        selected = selections[model_id]["fibrosis_group3"]
        if (
            selected.get("fallback_triggered") is not True
            or selected.get("valid_candidate_count") != 0
            or selected.get("candidate_failure_count") != 12
            or selected.get("outer_test_outcomes_used") is not False
            or selected.get("outer_test_features_used") is not False
            or selected.get("randomness_used") is not False
        ):
            raise SystemExit(f"v3 fallback selection differs: {path}: {model_id}")
    observed[(outer, seed)] = {
        "activation_count": activation_count,
        "activated_model_ids": activated,
    }
if set(observed) != set(expected):
    raise SystemExit("v3 validation unit roster differs")
(root / "validation/receipt.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-gse267145-production-v3-validation-v1",
            "status": "passed_unscored",
            "revision_id": contract["revision_id"],
            "unit_tests_passed": True,
            "v2_outer_1_seed_1701_predictions_byte_identical": True,
            "units": {f"{outer}/{seed}": value for (outer, seed), value in observed.items()},
            "primary_endpoint_code_path_changed": False,
            "other_secondary_endpoint_code_paths_changed": False,
            "outer_test_outcomes_read": False,
            "metrics_calculated": False,
            "scorer_called": False,
            "production_submission_authorized": False,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
PY

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  config/evaluation/gse267145_histology_production_v3.json \
  scripts/fit_gse267145_histology_baselines.py \
  scripts/fit_gse267145_histology_baselines_v3.py \
  tests/unit/test_fit_gse267145_histology_baselines.py \
  tests/unit/test_fit_gse267145_histology_baselines_v3.py \
  slurm/validate_gse267145_histology_production_v3_cpu.sbatch \
  "${MOLECULAR}/ARTIFACTS.json" \
  "${FOLDS}/ARTIFACTS.json" \
  "${FIT_VIEWS}/ARTIFACTS.json" \
  "${SURFACE}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${V2}/campaign_spec/ARTIFACTS.json" \
  "${V2}/units/outer_1/seed_1721.failed-21083268/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/receipt.json").read_text(encoding="utf-8"))
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_production_v3_validation",
        "logical_units": 2,
        "fallback_activations": sum(
            value["activation_count"] for value in receipt["units"].values()
        ),
        "outcomes_read": False,
        "metrics_calculated": False,
        "scorer_called": False,
        "production_submission_authorized": False,
        "status": "passed_unscored",
    },
)
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
