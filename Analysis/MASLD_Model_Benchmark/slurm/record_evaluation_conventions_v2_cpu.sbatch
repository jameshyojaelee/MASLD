#!/usr/bin/env bash
#SBATCH --job-name=model-check-904
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
DIAG=${ROOT}/executions/model-check-902-21177573-gain-concentration
DIAG_SHA256=c088a98fe1f3543f7d81a95189eb9fe412ac2e53e192de1773342e494d06c2d3
TARGET=${ROOT}/executions/model-check-904-${SLURM_JOB_ID}-evaluation-conventions
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite convention record: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}/convention" "${STAGE}/validation"
cd "${ROOT}"
hostname; "${PYTHON}" -V

printf '%s  %s\n' "${DIAG_SHA256}" "${DIAG}/ARTIFACTS.json" | sha256sum --check --strict

# Snapshot everything this record reads, so "modifies nothing" is evidence.
sha256sum \
  "${DIAG}/ARTIFACTS.json" \
  config/tasks/typed_evidence_graph.toml \
  config/evaluation/bulk_lineage_composition_fibrosis_promotion_gate.json \
  src/masld_bench/evaluators/auprc_reference.py \
  src/masld_bench/evaluators/stats.py \
  > "${STAGE}/validation/inputs_before.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_auprc_reference tests.unit.test_stats \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

# Regenerate the demonstration inside the job so it is reproducible from the
# frozen source rather than pasted from a session.
PYTHONPATH="${ROOT}/src" "${PYTHON}" - > "${STAGE}/convention/demonstration.json" <<'PY'
import json
from masld_bench.evaluators.auprc_reference import (
    permutation_reference, stratified_permutation_reference)
labels, scores, strata = [], [], []
for s, (rate, off) in enumerate(((0.9, 10.0), (0.7, 6.0), (0.3, 3.0), (0.1, 0.0))):
    for i in range(20):
        labels.append(1 if i < int(round(rate * 20)) else 0)
        scores.append(off + (i % 5) * 0.01)
        strata.append(f"s{s}")
g, obs, gp = permutation_reference(labels, scores, n_permutations=10000, seed=11)
w, _, wp = stratified_permutation_reference(labels, scores, strata, n_permutations=10000, seed=11)
c, _, _ = stratified_permutation_reference(
    labels, scores, [f"r{i % 4}" for i in range(len(labels))],
    n_permutations=10000, seed=11)
print(json.dumps({
    "fixture": ("four strata whose membership predicts the label; within each "
                "stratum the score carries NO information about the label, so "
                "an honest null must report no signal"),
    "n": len(labels), "n_strata": 4, "n_permutations": 10000, "seed": 11,
    "observed_average_precision": obs,
    "global_null_mean": g.mean, "global_null_p95": g.percentile_95,
    "global_p_value": gp,
    "global_p_value_note": "at the 1/(n+1) resolution floor; report as <= 1e-4",
    "stratified_null_mean": w.mean, "stratified_null_p95": w.percentile_95,
    "stratified_p_value": wp,
    "control_random_strata_null_mean": c.mean,
    "reading": ("Same data, same observed value. The global permutation "
                "destroys the between-stratum structure along with the signal "
                "and reports a highly significant result for a model that "
                "learned nothing within stratum. The stratified null correctly "
                "reports nothing. The control confirms the gap comes from real "
                "structure: with strata assigned at random the two nulls "
                "converge."),
}, indent=2, sort_keys=True))
PY

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" scripts/record_evaluation_conventions_v2.py \
  --demonstration "${STAGE}/convention/demonstration.json" \
  --retrospective "${DIAG}/diagnostic/gain_concentration.json" \
  --output "${STAGE}/convention/evaluation_conventions.json" \
  > "${STAGE}/validation/record_stdout.json"

sha256sum \
  "${DIAG}/ARTIFACTS.json" \
  config/tasks/typed_evidence_graph.toml \
  config/evaluation/bulk_lineage_composition_fibrosis_promotion_gate.json \
  src/masld_bench/evaluators/auprc_reference.py \
  src/masld_bench/evaluators/stats.py \
  > "${STAGE}/validation/inputs_after.sha256"
diff "${STAGE}/validation/inputs_before.sha256" "${STAGE}/validation/inputs_after.sha256"
printf 'every read input byte-identical before and after\n'

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum \
  scripts/record_evaluation_conventions_v2.py \
  src/masld_bench/evaluators/auprc_reference.py \
  src/masld_bench/evaluators/stats.py \
  slurm/record_evaluation_conventions_v2_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
record = json.loads((root / "convention" / "evaluation_conventions.json").read_text())
assert record["modifies_nothing"] is True
assert record["record_is_additive_overlay"] is True
assert "FORWARD-ONLY" in record["scope"]
assert record["convention_1_stratified_null"]["force"] == "PROCEDURAL REQUIREMENT"
assert "never pass/fail" in record["convention_2_gain_concentration"]["force"]
demo = record["convention_1_stratified_null"]["demonstration"]
assert demo["stratified_p_value"] > demo["global_p_value"], "demonstration must show the gap"
retro = record["convention_2_gain_concentration"]["retrospective_finding"]
assert retro["is_post_hoc"] is True and retro["revises_nothing"] is True
freeze_tree(root, {
    "artifact_class": "campaign_evaluation_convention",
    "record_is_additive_overlay": True,
    "modifies_nothing": True,
    "enforcement": "forward_only",
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'convention: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/convention/ARTIFACTS.json" 2>/dev/null || sha256sum "${TARGET}/ARTIFACTS.json"
