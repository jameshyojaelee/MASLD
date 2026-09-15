#!/usr/bin/env bash
#SBATCH --job-name=model-check-076
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/model-check-076-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
DIAGNOSTIC=${ROOT}/executions/model-check-075-21088989

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 v4 bundle validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

for pair in \
  "6ef4b1c6c460568fb982bfb174dd6047b69aad1a50931623b21b89c04df35d32 scripts/run_gse267145_histology_production_bundle_v4.py" \
  "1f4faa747eb936b61be1852f0279bcabb7021474f9bddbc7622e5628785889c9 slurm/run_gse267145_histology_production_bundle_v4_cpu.sbatch" \
  "b919d36182e484c0b8ae7b2f82043259c9e192cfca085713136ec945684821f7 tests/unit/test_gse267145_histology_production_bundle_v4_contract.py" \
  "df1f77f82c3657e6786be9a67100b3ce138384f8b28eb6a22b48bc3e0cf97de2 config/evaluation/gse267145_histology_production_v4.json" \
  "1e0f9ec881a3d04601f0e2c4c375458ef7664c72cd5fb733e400ab08e51f3cff scripts/fit_gse267145_histology_baselines_v4.py" \
  "0fd043f04f85328a3c06f67d2c65ae7f1115195e1c7a47be21a1b8708d21502a tests/unit/test_fit_gse267145_histology_baselines_v4.py" \
  "8616379cdf1f892fd994b754758028c38bf3d74f5b4dc20657848d3b7a8aedef config/evaluation/gse267145_histology_production_v3_reuse_manifest.json" \
  "b068e3d9af7952ce6bfe49d7482b4453c9bbb15e6d49cbfd3c912dba9b4fd60d executions/model-check-075-21088989/ARTIFACTS.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

bash -n slurm/run_gse267145_histology_production_bundle_v4_cpu.sbatch
"${PYTHON}" -m py_compile \
  scripts/run_gse267145_histology_production_bundle.py \
  scripts/run_gse267145_histology_production_bundle_v3.py \
  scripts/run_gse267145_histology_production_bundle_v4.py \
  scripts/fit_gse267145_histology_baselines_v4.py
env \
  PYTHONPATH="${ROOT}/src" \
  PYTHONDONTWRITEBYTECODE=1 \
  PYTHONNOUSERSITE=1 \
  OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  "${PYTHON}" scripts/fit_gse267145_histology_baselines_v4.py \
    --software-self-test \
    > "${STAGE}/validation/actual_child_environment_self_test.json"
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse267145_histology_production_bundle_contract \
  tests.unit.test_gse267145_histology_production_bundle_v3_contract \
  tests.unit.test_gse267145_histology_production_bundle_v4_contract \
  tests.unit.test_fit_gse267145_histology_baselines_v4 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" - "${STAGE}" <<'PY'
from argparse import Namespace
import json
import os
from pathlib import Path
import sys

import scripts.run_gse267145_histology_production_bundle_v4 as runner
from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
root = Path(runner.__file__).resolve().parents[1]
runtime = root / "executions/environments/transcriptformer-torch2.5.1-21070738"
preflight = root / "executions/model-training-067-21080776"
runner.v3.v2.FITTER_SHA256 = runner.V4_FITTER_SHA256
runner.v3.v2.DIAGNOSTIC_ARTIFACTS_SHA256 = runner.DIAGNOSTIC_ARTIFACTS_SHA256
self_test = json.loads(
    (stage / "validation/actual_child_environment_self_test.json").read_text(
        encoding="utf-8"
    )
)
if (
    self_test.get("status") != "passed_child_import_contract"
    or self_test.get("software_revision_id") != runner.REVISION_ID
    or self_test.get("repository_root_on_sys_path") is not True
    or self_test.get("v3_fallback_contract_sha256") != runner.v3.CONTRACT_SHA256
    or self_test.get("outcomes_read") is not False
    or self_test.get("metrics_calculated") is not False
    or self_test.get("scorer_called") is not False
):
    raise SystemExit("v4 child-environment self-test differs")

production_receipt = {
    "schema_version": "masld-bench-gse267145-production-v4-bundle-validation-v1",
    "status": "passed_unscored_production_wrapper_not_submitted",
    "revision_id": runner.REVISION_ID,
    "parity_validation_artifacts_sha256": runner.DIAGNOSTIC_ARTIFACTS_SHA256,
    "reused_v2_units_byte_identical": 22,
    "v4_units_to_compute": ["1/1721", "1/1723", "1/1733"],
    "failed_v3_attempts_reused": False,
    "participant_outer_folds_artifacts_sha256": runner.v3.v2.FOLDS_ARTIFACTS_SHA256,
    "production_slurm": {
        "job_name": "model-training-077",
        "partition": "cpu",
        "account": "nslab",
        "qos": "nslab",
        "cpus": 16,
        "memory": "32G",
        "wall_time": "04:00:00",
        "fold_workers": 5,
        "blas_threads_per_active_fit": 3,
        "maximum_concurrent_new_unit_fits_after_reuse": 1,
        "maximum_cpu_hours": 64,
        "maximum_memory_gb_hours": 128,
    },
    "actual_child_environment_self_test_passed": True,
    "outer_test_outcomes_read": False,
    "metrics_calculated": False,
    "scorer_called": False,
    "production_wrapper_submitted": False,
}
binding_fixture = stage / "binding_fixture"
(binding_fixture / "validation").mkdir(parents=True)
(binding_fixture / "validation/receipt.json").write_text(
    json.dumps(production_receipt, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
required_sources = (
    "scripts/run_gse267145_histology_production_bundle.py",
    "scripts/run_gse267145_histology_production_bundle_v3.py",
    "scripts/run_gse267145_histology_production_bundle_v4.py",
    "scripts/fit_gse267145_histology_baselines.py",
    "scripts/fit_gse267145_histology_baselines_v3.py",
    "scripts/fit_gse267145_histology_baselines_v4.py",
    "scripts/aggregate_gse267145_histology_production_predictions.py",
    "config/evaluation/gse267145_histology_production_v3.json",
    "config/evaluation/gse267145_histology_production_v4.json",
    "config/evaluation/gse267145_histology_production_v3_reuse_manifest.json",
    "slurm/run_gse267145_histology_production_bundle_v4_cpu.sbatch",
)
source_rows = [
    f"{runner.v3.v2.sha256_file(root / relative)}  {relative}"
    for relative in required_sources
]
diagnostic_artifacts = root / "executions/model-check-075-21088989/ARTIFACTS.json"
source_rows.append(
    f"{runner.v3.v2.sha256_file(diagnostic_artifacts)}  {diagnostic_artifacts.as_posix()}"
)
(binding_fixture / "source.sha256").write_text(
    "\n".join(source_rows) + "\n", encoding="utf-8"
)
binding_sha256 = freeze_tree(
    binding_fixture,
    {
        "artifact_class": "gse267145_histology_production_v4_bundle_validation",
        "reused_v2_units_byte_identical": 22,
        "v4_units_to_compute": 3,
        "outcomes_read": False,
        "metrics_calculated": False,
        "scorer_called": False,
        "production_wrapper_submitted": False,
        "status": "passed_unscored",
    },
)
os.environ["BUNDLE_VALIDATION_ROOT"] = binding_fixture.as_posix()
os.environ["BUNDLE_VALIDATION_SHA256"] = binding_sha256
arguments = Namespace(
    runtime=runtime,
    molecular=root / "executions/model-data-064-21079902/fixture/molecular",
    folds=root / "executions/model-data-064-21079902/fixture/folds",
    fit_views=preflight / "fit_views",
    preflight=preflight,
    diagnostic=root / "executions/model-check-075-21088989",
    surface=root / "config/evaluation/gse267145_histology_benchmark_surface.json",
    fitter=root / "scripts/fit_gse267145_histology_baselines_v4.py",
    python=runtime / "env/bin/python",
    fold_workers=5,
    blas_threads=3,
    campaign=root / "executions/model-training-069-production-v4",
)
runner._input_contract_v4(arguments)

fixture = stage / "reuse_fixture"
runner.v3._populate_reused_units(Namespace(campaign=fixture), "validation-v4")
reused, recomputed = runner.v3._reuse_manifest()
observed = set()
for path in sorted((fixture / "units").glob("outer_*/seed_*")):
    outer = int(path.parent.name.removeprefix("outer_"))
    seed = int(path.name.removeprefix("seed_"))
    key = (outer, seed)
    runner.v3._verify_reused_unit(path, outer, seed, reused[key])
    observed.add(key)
if observed != set(reused) or len(observed) != 22:
    raise SystemExit("v4 copied reuse roster differs")
if recomputed != {(1, 1721), (1, 1723), (1, 1733)}:
    raise SystemExit("v4 recomputation roster differs")
if any(
    (fixture / "units" / f"outer_{outer}" / f"seed_{seed}").exists()
    for outer, seed in recomputed
):
    raise SystemExit("v4 reuse fixture populated a recomputation unit")

(stage / "validation/receipt.json").write_text(
    json.dumps(production_receipt, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/run_gse267145_histology_production_bundle.py \
  scripts/run_gse267145_histology_production_bundle_v3.py \
  scripts/run_gse267145_histology_production_bundle_v4.py \
  scripts/fit_gse267145_histology_baselines.py \
  scripts/fit_gse267145_histology_baselines_v3.py \
  scripts/fit_gse267145_histology_baselines_v4.py \
  scripts/aggregate_gse267145_histology_production_predictions.py \
  config/evaluation/gse267145_histology_production_v3.json \
  config/evaluation/gse267145_histology_production_v4.json \
  config/evaluation/gse267145_histology_production_v3_reuse_manifest.json \
  tests/unit/test_gse267145_histology_production_bundle_contract.py \
  tests/unit/test_gse267145_histology_production_bundle_v3_contract.py \
  tests/unit/test_gse267145_histology_production_bundle_v4_contract.py \
  tests/unit/test_fit_gse267145_histology_baselines_v4.py \
  slurm/run_gse267145_histology_production_bundle_v4_cpu.sbatch \
  slurm/validate_gse267145_histology_production_bundle_v4_cpu.sbatch \
  "${DIAGNOSTIC}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_production_v4_bundle_validation",
        "reused_v2_units_byte_identical": 22,
        "v4_units_to_compute": 3,
        "actual_child_environment_self_test_passed": True,
        "outcomes_read": False,
        "metrics_calculated": False,
        "scorer_called": False,
        "production_wrapper_submitted": False,
        "status": "passed_unscored",
    },
)
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
