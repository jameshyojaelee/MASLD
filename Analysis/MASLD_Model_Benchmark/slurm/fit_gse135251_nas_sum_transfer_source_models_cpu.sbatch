#!/usr/bin/env bash
#SBATCH --job-name=model-fit-931
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
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
RUNTIME_SHA=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CONTRACT=${ROOT}/config/evaluation/gse135251_nas_sum_transfer_source_fit_v1.json
CONTRACT_SHA=06a5e18520c85ba1afd134ea69db1a275623384289447ab847bfa50304c7925b
ACTIVATION=${ROOT}/executions/model-data-929-21133561/activation
SOURCE=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
OUT=${ROOT}/executions/model-fit-931-${SLURM_JOB_ID}
STAGE=${OUT}.staging
FAILED=${OUT}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${OUT}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite NAS source-fit artifacts.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${CONTRACT_SHA}" "${CONTRACT}" | sha256sum --check --strict
bash -n slurm/fit_gse135251_nas_sum_transfer_source_models_cpu.sbatch
# The source fit may not name the target outcome table or its histology columns.
if grep -n 'participant_endpoints\.tsv.*267145\|nash_crn_component_sum\|model-data-061' \
    scripts/fit_gse135251_nas_sum_transfer_source_models.py; then
  printf 'Source fit unexpectedly names the target outcome.\n' >&2; exit 1
fi
"${PYTHON}" -m unittest -v \
  tests.unit.test_fit_gse135251_nas_sum_transfer_source_models \
  tests.unit.test_fit_gse267145_fibrosis_transfer_source_models \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/fit_gse135251_nas_sum_transfer_source_models.py \
  --benchmark-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --contract-sha256 "${CONTRACT_SHA}" \
  --output "${STAGE}/source_fit" \
  > "${STAGE}/fit_stdout.json" 2> "${STAGE}/fit_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/fit_gse135251_nas_sum_transfer_source_models.py \
  scripts/fit_gse267145_fibrosis_transfer_source_models.py \
  tests/unit/test_fit_gse135251_nas_sum_transfer_source_models.py \
  slurm/fit_gse135251_nas_sum_transfer_source_models_cpu.sbatch \
  "${CONTRACT}" "${RUNTIME}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" "${SOURCE}/ARTIFACTS.json" \
  config/evaluation/gse267145_nas_sum_transfer_task.toml \
  config/evaluation/gse267145_nas_sum_transfer_promotion_gate.json \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv, json, sys
import numpy as np
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
fit = root / "source_fit"
r = json.loads((fit / "fit_receipt.json").read_text(encoding="utf-8"))
with (fit / "seed_model_index.tsv").open(encoding="utf-8", newline="") as h:
    seeds = list(csv.DictReader(h, delimiter="\t"))
if (
    r.get("status") != "passed_source_nas_sum_fit_preprocessing_locked"
    or r.get("source_participants") != 180
    or r.get("shared_genes") != 38787
    or r.get("endpoint_binarised") is not False
    or r.get("target_outcomes_read") is not False
    or r.get("target_expression_values_read") is not False
    or r.get("cross_sample_pooling_performed") is not False
    or len(seeds) != 5 * len(r["fitted_model_ids"])
    or not r["fitted_model_ids"]
):
    raise SystemExit("NAS source-fit verification differs")
for row in seeds:
    if int(row["nonzero_coefficients"]) == 0:
        raise SystemExit(f"degenerate model reached the freeze: {row['model_id']}")
    seed_root = (fit / row["model_state_path"]).parent
    freeze_tree(seed_root, {
        "artifact_class": "gse135251_nas_sum_seed_model",
        "model_id": row["model_id"], "model_seed": int(row["model_seed"]),
        "target_data_read": False, "status": "passed",
    })
    verify_frozen_tree(seed_root)
freeze_tree(fit, {
    "artifact_class": "gse135251_nas_sum_source_fit",
    "source_participants": 180, "shared_genes": 38787,
    "fitted_models": len(r["fitted_model_ids"]),
    "selected_source_model_id": r["selected_source_model_id"],
    "target_outcomes_read": False, "status": "passed",
})
verify_frozen_tree(fit)
freeze_tree(root, {
    "artifact_class": "gse135251_nas_sum_source_fit_campaign",
    "target_outcomes_read": False, "metrics_calculated": False, "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "selected_source_model_id": r["selected_source_model_id"],
    "fitted_model_ids": r["fitted_model_ids"],
    "unfittable_model_ids": r["unfittable_model_ids"],
    "no_source_model_exceeds_its_own_null": r["no_source_model_exceeds_its_own_null"],
    "training_mean_nas_sum": r["training_mean_nas_sum"],
    "source_model_ranking": r["source_model_ranking"],
}, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${OUT}"
sha256sum "${OUT}/ARTIFACTS.json" "${OUT}/source_fit/ARTIFACTS.json"
