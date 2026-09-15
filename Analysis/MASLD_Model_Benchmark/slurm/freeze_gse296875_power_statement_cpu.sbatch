#!/usr/bin/env bash
#SBATCH --job-name=model-check-011
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
TARGET="${ROOT}/executions/gse296875-power-statement-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite power-statement artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
mkdir -p "${STAGE}/validation"

cat > "${STAGE}/validation/bound_fixtures.sha256" <<'SUMS'
6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140  executions/gse296875-phenotype-endpoints-20260825/ARTIFACTS.json
e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416  executions/gse296875-donor-folds-20260825/ARTIFACTS.json
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
3f371ba28847f28c5c17295d50dc6ca65f51c4bde6039ee31e2bcf764603d40f  executions/ranking-metric-null-reference-note-v3-21109019/ARTIFACTS.json
d1343044e312dd88c1026f4f6af1b6f284d20f4f4b8e155c6fb3e45964d956d8  executions/gse296875-decision-thresholds-v2-21114548/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

PYTHONPATH=src python3 -m unittest -v \
  tests.unit.test_detectable_effect \
  tests.unit.test_auprc_reference \
  tests.contract.test_gse296875_phenotype \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH=src python3 scripts/freeze_gse296875_power_statement.py \
  --root "${ROOT}" \
  --output "${STAGE}/power" \
  --realised \
    "${ROOT}/executions/gse296875-realised-power-steatosis-vf0.02-21113440" \
    "${ROOT}/executions/gse296875-realised-power-steatosis-vf0.10-21113441" \
    "${ROOT}/executions/gse296875-realised-power-fibrosis-vf0.02-21113442" \
    "${ROOT}/executions/gse296875-realised-power-fibrosis-vf0.10-21113443" \
  > "${STAGE}/validation/statement_stdout.json"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum \
  scripts/freeze_gse296875_power_statement.py \
  src/masld_bench/evaluators/detectable_effect.py \
  tests/unit/test_detectable_effect.py \
  slurm/freeze_gse296875_power_statement_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
statement = json.loads((root / "power" / "power_statement.json").read_text())
if statement["observed_model_performance_read"] is not False:
    raise SystemExit("a power statement may not read an observed result")
if statement["derived_from"] != "frozen_design_only":
    raise SystemExit("the statement must be derived from the design alone")
design = statement["design"]
if design["steatosis_observed"] != 38 or design["fibrosis_observed"] != 37:
    raise SystemExit("frozen denominators differ")
if design["fibrosis_positive"] != 15 or design["fibrosis_negative"] != 22:
    raise SystemExit("frozen class balance differs")
if design["fold_3_has_zero_fibrosis_positives"] is not True:
    raise SystemExit("fold 3 fact missing")
if design["pooled_out_of_fold_metric_only"] is not True:
    raise SystemExit("pooled-only requirement missing")
freeze_tree(root, {
    "artifact_class": "gse296875_prespecified_power_statement",
    "dataset_id": "gse296875",
    "unit_of_inference": "donor",
    "derived_from": "frozen_design_only",
    "observed_model_performance_read": False,
    "campaign_spec_sha256": "474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58",
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE296875 power statement: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
