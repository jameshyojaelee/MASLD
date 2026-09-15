#!/usr/bin/env bash
#SBATCH --job-name=model-check-1043
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 0 of the MASLD showcase model, job 1 of 2.
#
# Freezes what would count as three separable disease aspects BEFORE the
# correlations that answer the question are computed. The analysis job
# (model-check-1044) verifies this tree's ARTIFACTS.json digest before it
# computes anything, so the ordering is provable from the output files alone and
# not from a claim in a report.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/freeze_aspect_separability_prespecification.py
ENDPOINTS=${ROOT}/executions/model-data-061-21079623/activation/participant_endpoints.tsv
ENDPOINTS_SHA256=9a13c0fd65db0483566b7495def11ade3339861de091ea29fcbdc4c34f47b5b9
TARGET=${ROOT}/executions/model-check-1043-${SLURM_JOB_ID}-aspect-separability-prespec
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the prespecification record.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/freeze_aspect_separability_prespecification_cpu.sbatch

# The endpoint table is the only frozen input this job reads. Its marginals are
# what the prespecification records as inspected before the freeze, so the
# digest is checked before and after.
printf '%s  %s\n' "${ENDPOINTS_SHA256}" "${ENDPOINTS}" | sha256sum --check --strict

"${PYTHON}" -m unittest -v tests.unit.test_freeze_aspect_separability_prespecification \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

# The script owns its --output path and refuses to overwrite one. The wrapper
# must not pre-create the directory it is handed.
"${PYTHON}" "${SCRIPT}" \
  --output "${STAGE}/prespecification/aspect_separability_prespecification.json" \
  > "${STAGE}/prespecification_stdout.json" 2> "${STAGE}/prespecification_stderr.txt"

printf '%s  %s\n' "${ENDPOINTS_SHA256}" "${ENDPOINTS}" | sha256sum --check --strict

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" \
  slurm/freeze_aspect_separability_prespecification_cpu.sbatch \
  tests/unit/test_freeze_aspect_separability_prespecification.py \
  "${ENDPOINTS}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
prespec_dir = root / "prespecification"
payload = json.loads(
    (prespec_dir / "aspect_separability_prespecification.json").read_text(
        encoding="utf-8"
    )
)

expected = [
    "c1_labels_are_not_redundant",
    "c2_effective_dimensionality",
    "c3_aspects_rank_genes_differently",
]
if sorted(payload["criteria"]) != expected:
    raise SystemExit("the three criteria are not the frozen three")
if payload["criteria"]["c1_labels_are_not_redundant"]["threshold"] != (
    "every pair |rho| < 0.80"
):
    raise SystemExit("c1 threshold differs")
if payload["criteria"]["c2_effective_dimensionality"]["threshold"] != (
    ">= 2.0 of a possible 3.0"
):
    raise SystemExit("c2 threshold differs")
if payload["criteria"]["c3_aspects_rank_genes_differently"]["threshold"] != (
    "at least one pair with |rho| < 0.80"
):
    raise SystemExit("c3 threshold differs")
if payload["criteria"]["c3_aspects_rank_genes_differently"].get("decisive") is not True:
    raise SystemExit("c3 is not marked decisive")
if payload["decision_rule"]["a_stop_is_a_result_not_a_failure"] is not True:
    raise SystemExit("a stop is not recorded as a result")
inspected = payload["what_was_inspected_before_freezing"]
if inspected["no_gene_expression_was_opened"] is not True:
    raise SystemExit("the prespecification does not assert an unopened matrix")
if inspected["nash_crn_component_sum_zero_participants"] != 24:
    raise SystemExit("the recorded all-zero block differs")

freeze_tree(prespec_dir, {
    "artifact_class": "gse267145_aspect_separability_prespecification",
    "prespec_id": payload["prespec_id"],
    "criteria_count": 3,
    "decisive_criterion": "c3_aspects_rank_genes_differently",
    "computed_before_any_aspect_correlation": True,
    "gene_expression_opened": False,
    "status": "passed",
})
verify_frozen_tree(prespec_dir)
freeze_tree(root, {
    "artifact_class": "gse267145_aspect_separability_prespec_campaign",
    "stage": "stage_0_gate",
    "job": "1_of_2",
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "prespec_id": payload["prespec_id"],
    "criteria": sorted(payload["criteria"]),
    "decision_rule_go": payload["decision_rule"]["go"],
}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/prespecification/ARTIFACTS.json"
