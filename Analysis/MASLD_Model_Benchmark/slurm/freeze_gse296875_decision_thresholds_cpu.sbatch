#!/usr/bin/env bash
#SBATCH --job-name=model-check-014
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/gse296875-decision-thresholds-v2-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite threshold artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
mkdir -p "${STAGE}/validation"

cat > "${STAGE}/validation/bound_fixtures.sha256" <<'SUMS'
6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140  executions/gse296875-phenotype-endpoints-20260825/ARTIFACTS.json
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

PYTHONPATH=src python3 -m unittest -v \
  tests.unit.test_gse296875_phenotype_scoring \
  tests.unit.test_auprc_reference \
  tests.unit.test_detectable_effect \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH=src python3 scripts/freeze_gse296875_decision_thresholds.py \
  --root "${ROOT}" --output "${STAGE}/thresholds" \
  > "${STAGE}/validation/threshold_stdout.json"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum \
  scripts/freeze_gse296875_decision_thresholds.py \
  src/masld_bench/evaluators/detectable_effect.py \
  src/masld_bench/gse296875_phenotype_scoring.py \
  slurm/freeze_gse296875_decision_thresholds_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "thresholds" / "decision_thresholds.json").read_text())
if receipt["observed_model_performance_read"] is not False:
    raise SystemExit("thresholds may not be derived from an observed result")
if receipt["derived_inside_the_frozen_tree"] is not True:
    raise SystemExit("thresholds must be derived inside the tree")
for key, record in receipt["thresholds"].items():
    if record["total_draws"] < 500_000:
        raise SystemExit(f"{key} needs at least 500,000 draws")
    if record["percentile_scale"] not in {"absolute_two_sided", "signed_one_sided"}:
        raise SystemExit(f"{key} does not declare its side")
    for field in ("label_vector_ordering", "quantile_convention", "seed"):
        if not record.get(field):
            raise SystemExit(f"{key} does not record {field}")
if "label_vector_ordering_is_a_free_parameter" not in receipt:
    raise SystemExit("the ordering parameter must be documented")
if not receipt["label_vector_ordering_is_a_free_parameter"]["demonstration"]:
    raise SystemExit("the ordering effect must be demonstrated, not asserted")
freeze_tree(root, {
    "artifact_class": "gse296875_decision_thresholds",
    "dataset_id": "gse296875",
    "derived_inside_the_frozen_tree": True,
    "total_draws_per_threshold": 1_000_000,
    "records_all_four_estimate_parameters": True,
    "observed_model_performance_read": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE296875 decision thresholds: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
