#!/usr/bin/env bash
#SBATCH --job-name=model-check-074
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
TARGET=${ROOT}/executions/model-check-074-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
DIAGNOSTIC=${ROOT}/executions/model-check-072-21088262

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 v3 bundle validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

for pair in \
  "1fd23e06c0ea5b18d1c0e0a76ba82bc75b5237c7a1f63b149552381f17aee54c scripts/run_gse267145_histology_production_bundle_v3.py" \
  "6829afd4d118c7286b483f6fdcda8e7c8e77e7faf1da080c09a973e70512465d slurm/run_gse267145_histology_production_bundle_v3_cpu.sbatch" \
  "76ac05f1fbe3876afad64d887bc765f38b9361573632fcf979d3295c575f2df6 tests/unit/test_gse267145_histology_production_bundle_v3_contract.py" \
  "8616379cdf1f892fd994b754758028c38bf3d74f5b4dc20657848d3b7a8aedef config/evaluation/gse267145_histology_production_v3_reuse_manifest.json" \
  "8be11a742def0f80465bca92fd02baf626d2d207ee30a87cdb25bcf59e47e29f config/evaluation/gse267145_histology_production_v3.json" \
  "e65a54a3afefd386ae4fa239de37665ad45bf7b8ce34f7ad8f0a15157653972a executions/model-check-072-21088262/ARTIFACTS.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

bash -n slurm/run_gse267145_histology_production_bundle_v3_cpu.sbatch
"${PYTHON}" -m py_compile \
  scripts/run_gse267145_histology_production_bundle.py \
  scripts/run_gse267145_histology_production_bundle_v3.py
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse267145_histology_production_bundle_contract \
  tests.unit.test_gse267145_histology_production_bundle_v3_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" - "${STAGE}" <<'PY'
from argparse import Namespace
import json
import os
from pathlib import Path
import sys

import scripts.run_gse267145_histology_production_bundle_v3 as runner
from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
root = Path(runner.__file__).resolve().parents[1]
runtime = root / "executions/environments/transcriptformer-torch2.5.1-21070738"
preflight = root / "executions/model-training-067-21080776"
runner.v2.FITTER_SHA256 = runner.V3_FITTER_SHA256
runner.v2.DIAGNOSTIC_ARTIFACTS_SHA256 = runner.DIAGNOSTIC_ARTIFACTS_SHA256
production_receipt = {
    "schema_version": "masld-bench-gse267145-production-v3-bundle-validation-v1",
    "status": "passed_unscored_production_wrapper_not_submitted",
    "revision_id": runner.REVISION_ID,
    "fallback_validation_artifacts_sha256": runner.DIAGNOSTIC_ARTIFACTS_SHA256,
    "reused_v2_units_byte_identical": 22,
    "v3_units_to_compute": ["1/1721", "1/1723", "1/1733"],
    "failed_v2_attempt_reused": False,
    "participant_outer_folds_artifacts_sha256": runner.v2.FOLDS_ARTIFACTS_SHA256,
    "production_slurm": {
        "job_name": "model-training-073",
        "partition": "cpu",
        "cpus": 16,
        "memory": "32G",
        "wall_time": "04:00:00",
        "fold_workers": 5,
        "blas_threads_per_active_fit": 3,
        "maximum_concurrent_new_unit_fits_after_reuse": 1,
    },
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
    "scripts/fit_gse267145_histology_baselines.py",
    "scripts/fit_gse267145_histology_baselines_v3.py",
    "scripts/aggregate_gse267145_histology_production_predictions.py",
    "config/evaluation/gse267145_histology_production_v3.json",
    "config/evaluation/gse267145_histology_production_v3_reuse_manifest.json",
    "slurm/run_gse267145_histology_production_bundle_v3_cpu.sbatch",
)
source_rows = [
    f"{runner.v2.sha256_file(root / relative)}  {relative}"
    for relative in required_sources
]
diagnostic_artifacts = root / "executions/model-check-072-21088262/ARTIFACTS.json"
source_rows.append(
    f"{runner.v2.sha256_file(diagnostic_artifacts)}  {diagnostic_artifacts.as_posix()}"
)
(binding_fixture / "source.sha256").write_text(
    "\n".join(source_rows) + "\n", encoding="utf-8"
)
binding_sha256 = freeze_tree(
    binding_fixture,
    {
        "artifact_class": "gse267145_histology_production_v3_bundle_validation",
        "reused_v2_units_byte_identical": 22,
        "v3_units_to_compute": 3,
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
    diagnostic=root / "executions/model-check-072-21088262",
    surface=root / "config/evaluation/gse267145_histology_benchmark_surface.json",
    fitter=root / "scripts/fit_gse267145_histology_baselines_v3.py",
    python=runtime / "env/bin/python",
    fold_workers=5,
    blas_threads=3,
    campaign=root / "executions/model-training-069-production-v3",
)
runner._input_contract_v3(arguments)

fixture = stage / "reuse_fixture"
runner._populate_reused_units(Namespace(campaign=fixture), "validation")
reused, recomputed = runner._reuse_manifest()
observed = set()
for path in sorted((fixture / "units").glob("outer_*/seed_*")):
    outer = int(path.parent.name.removeprefix("outer_"))
    seed = int(path.name.removeprefix("seed_"))
    key = (outer, seed)
    runner._verify_reused_unit(path, outer, seed, reused[key])
    observed.add(key)
if observed != set(reused) or len(observed) != 22:
    raise SystemExit("v3 copied reuse roster differs")
if recomputed != {(1, 1721), (1, 1723), (1, 1733)}:
    raise SystemExit("v3 recomputation roster differs")
if any((fixture / "units" / f"outer_{outer}" / f"seed_{seed}").exists()
       for outer, seed in recomputed):
    raise SystemExit("v3 reuse fixture populated a recomputation unit")

(stage / "validation/receipt.json").write_text(
    json.dumps(production_receipt, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/run_gse267145_histology_production_bundle.py \
  scripts/run_gse267145_histology_production_bundle_v3.py \
  scripts/fit_gse267145_histology_baselines.py \
  scripts/fit_gse267145_histology_baselines_v3.py \
  scripts/aggregate_gse267145_histology_production_predictions.py \
  config/evaluation/gse267145_histology_production_v3.json \
  config/evaluation/gse267145_histology_production_v3_reuse_manifest.json \
  tests/unit/test_gse267145_histology_production_bundle_contract.py \
  tests/unit/test_gse267145_histology_production_bundle_v3_contract.py \
  slurm/run_gse267145_histology_production_bundle_v3_cpu.sbatch \
  slurm/validate_gse267145_histology_production_bundle_v3_cpu.sbatch \
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
        "artifact_class": "gse267145_histology_production_v3_bundle_validation",
        "reused_v2_units_byte_identical": 22,
        "v3_units_to_compute": 3,
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
