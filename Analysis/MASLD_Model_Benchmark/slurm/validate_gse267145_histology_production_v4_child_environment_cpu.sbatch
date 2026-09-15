#!/usr/bin/env bash
#SBATCH --job-name=model-check-1042
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
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
CAMPAIGN=${ROOT}/executions/model-training-069-production-v4
V2=${ROOT}/executions/model-training-069-production-v2
PARITY=${ROOT}/executions/model-check-075-21088989
BUNDLE=${ROOT}/executions/model-check-076-21089609
TARGET=${ROOT}/executions/model-check-1042-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 v4 child-environment validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

for pair in \
  "df1f77f82c3657e6786be9a67100b3ce138384f8b28eb6a22b48bc3e0cf97de2 config/evaluation/gse267145_histology_production_v4.json" \
  "8be11a742def0f80465bca92fd02baf626d2d207ee30a87cdb25bcf59e47e29f config/evaluation/gse267145_histology_production_v3.json" \
  "8616379cdf1f892fd994b754758028c38bf3d74f5b4dc20657848d3b7a8aedef config/evaluation/gse267145_histology_production_v3_reuse_manifest.json" \
  "1e0f9ec881a3d04601f0e2c4c375458ef7664c72cd5fb733e400ab08e51f3cff scripts/fit_gse267145_histology_baselines_v4.py" \
  "0cf28366d86040ce18fec73e3a5cde862ab33dad71cf0719bdd1d15a14508d3e scripts/fit_gse267145_histology_baselines_v3.py" \
  "5e2caef9396a040173842df59d9c01e83cd6ac35031fe7050045e6209ad3ca8e scripts/run_gse267145_histology_production_bundle.py" \
  "1fd23e06c0ea5b18d1c0e0a76ba82bc75b5237c7a1f63b149552381f17aee54c scripts/run_gse267145_histology_production_bundle_v3.py" \
  "6ef4b1c6c460568fb982bfb174dd6047b69aad1a50931623b21b89c04df35d32 scripts/run_gse267145_histology_production_bundle_v4.py" \
  "b068e3d9af7952ce6bfe49d7482b4453c9bbb15e6d49cbfd3c912dba9b4fd60d executions/model-check-075-21088989/ARTIFACTS.json" \
  "5cac628ebe1714759288a50344829d605aafea37eb3a8d6042c7b5fda573849e executions/model-check-076-21089609/ARTIFACTS.json" \
  "ebbe02b9cf0e11ea6d01ccf981652631801b382fbd12225710b9fe17013da59d executions/model-training-069-production-v4/ARTIFACTS.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

"${PYTHON}" -m py_compile \
  tests/unit/test_gse267145_histology_production_child_environment.py
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse267145_histology_production_child_environment \
  tests.unit.test_fit_gse267145_histology_baselines_v4 \
  tests.unit.test_gse267145_histology_production_bundle_v4_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" - "${STAGE}" "${CAMPAIGN}" "${V2}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree
import scripts.run_gse267145_histology_production_bundle_v4 as runner

stage = Path(sys.argv[1]).resolve(strict=True)
campaign = Path(sys.argv[2]).resolve(strict=True)
v2_campaign = Path(sys.argv[3]).resolve(strict=True)

verify_frozen_tree(campaign)
verify_frozen_tree(campaign / "aggregate")
receipt = json.loads((campaign / "campaign_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("schema_version") != "masld-bench-gse267145-production-bundle-v1"
    or receipt.get("campaign_id") != runner.REVISION_ID
    or receipt.get("status") != "passed_predictions_unscored"
    or receipt.get("logical_outer_seed_units") != 25
    or receipt.get("model_family_unit_fits") != 275
    or receipt.get("outcomes_read") is not False
    or receipt.get("metrics_calculated") is not False
    or receipt.get("scorer_called") is not False
):
    raise SystemExit("v4 campaign receipt differs")

reused, recomputed = runner.v3._reuse_manifest()
if len(reused) != 22 or recomputed != {(1, 1721), (1, 1723), (1, 1733)}:
    raise SystemExit("v4 unit roster differs")

for (outer, seed), expected in sorted(reused.items()):
    unit = campaign / "units" / f"outer_{outer}" / f"seed_{seed}"
    runner.v3._verify_reused_unit(unit, outer, seed, expected)
    source = runner.v3._v2_unit(outer, seed)
    for prediction in sorted((source / "fit/predictions").glob("*.tsv")):
        target = unit / "fit/predictions" / prediction.name
        if prediction.read_bytes() != target.read_bytes():
            raise SystemExit(f"reused v2 unit is not byte-identical: {outer}/{seed}")

computed = {}
for outer, seed in sorted(recomputed):
    unit = campaign / "units" / f"outer_{outer}" / f"seed_{seed}"
    verify_frozen_tree(unit)
    unit_receipt = json.loads((unit / "unit_receipt.json").read_text(encoding="utf-8"))
    if (
        unit_receipt.get("schema_version")
        != "masld-bench-gse267145-production-unit-v1"
        or unit_receipt.get("status") != "passed_unscored"
        or unit_receipt.get("outer_fold") != outer
        or unit_receipt.get("seed") != seed
        or unit_receipt.get("production_unit_complete") is not True
        or unit_receipt.get("model_family_fits") != 11
        or unit_receipt.get("prediction_files") != 11
        or unit_receipt.get("outer_test_outcomes_read") is not False
        or unit_receipt.get("outer_test_metrics_calculated") is not False
        or unit_receipt.get("scorer_called") is not False
    ):
        raise SystemExit(f"recomputed v4 unit receipt differs: {outer}/{seed}")
    if (unit / "fit/predictions").exists():
        files = sorted(path.name for path in (unit / "fit/predictions").glob("*.tsv"))
        if len(files) != 11:
            raise SystemExit(f"recomputed v4 unit prediction count differs: {outer}/{seed}")
    computed[f"{outer}/{seed}"] = unit_receipt["prediction_artifacts_sha256"]

if any((v2_campaign / "units" / f"outer_{outer}" / f"seed_{seed}").exists()
       for outer, seed in recomputed):
    raise SystemExit("v2 campaign unexpectedly holds a recomputed unit")

(stage / "validation/receipt.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-gse267145-production-v4-child-environment-v1",
            "status": "passed_unscored",
            "revision_id": runner.REVISION_ID,
            "production_slurm_job_id": 21089724,
            "superseded_failed_slurm_jobs": [21088904, 21088907],
            "failure": "ModuleNotFoundError: No module named 'scripts'",
            "child_pythonpath": "<root>/src",
            "child_import_shim": "scripts/fit_gse267145_histology_baselines_v4.py",
            "campaign_artifacts_sha256": runner.v3.v2.sha256_file(
                campaign / "ARTIFACTS.json"
            ),
            "campaign_spec_artifacts_sha256": receipt["campaign_spec_artifacts_sha256"],
            "aggregate_artifacts_sha256": receipt["aggregate_artifacts_sha256"],
            "prediction_files": receipt["prediction_files"],
            "reused_v2_units_byte_identical": len(reused),
            "recomputed_v4_units": computed,
            "unit_tests_passed": True,
            "outer_test_outcomes_read": False,
            "metrics_calculated": False,
            "scorer_called": False,
            "evaluation_authorized": False,
            "scoring_authorized": False,
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
  config/evaluation/gse267145_histology_production_v4.json \
  config/evaluation/gse267145_histology_production_v3.json \
  config/evaluation/gse267145_histology_production_v3_reuse_manifest.json \
  scripts/fit_gse267145_histology_baselines.py \
  scripts/fit_gse267145_histology_baselines_v3.py \
  scripts/fit_gse267145_histology_baselines_v4.py \
  scripts/run_gse267145_histology_production_bundle.py \
  scripts/run_gse267145_histology_production_bundle_v3.py \
  scripts/run_gse267145_histology_production_bundle_v4.py \
  scripts/aggregate_gse267145_histology_production_predictions.py \
  tests/unit/test_gse267145_histology_production_child_environment.py \
  tests/unit/test_fit_gse267145_histology_baselines_v4.py \
  tests/unit/test_gse267145_histology_production_bundle_v4_contract.py \
  slurm/run_gse267145_histology_production_bundle_v4_cpu.sbatch \
  slurm/validate_gse267145_histology_production_v4_child_environment_cpu.sbatch \
  "${PARITY}/ARTIFACTS.json" \
  "${BUNDLE}/ARTIFACTS.json" \
  "${CAMPAIGN}/ARTIFACTS.json" \
  "${CAMPAIGN}/aggregate/ARTIFACTS.json" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_production_v4_child_environment",
        "reused_v2_units_byte_identical": 22,
        "recomputed_v4_units": 3,
        "child_import_regression_test_passed": True,
        "outcomes_read": False,
        "metrics_calculated": False,
        "scorer_called": False,
        "status": "passed_unscored",
    },
)
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
