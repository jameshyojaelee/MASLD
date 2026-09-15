#!/usr/bin/env bash
#SBATCH --job-name=model-check-899
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
W3REC=${ROOT}/executions/model-check-898-21162928-w3-recorded-negative/record
W3REC_SHA256=986e755992867bdb27a73af1845989d6e6e728270dd2a0380f7ed3ae67bd7626
W2GATE=${ROOT}/config/evaluation/bulk_lineage_composition_fibrosis_promotion_gate.json
W2GATE_SHA256=e9cb9f9d11d1023450b5a50115f0b1b61b326caf99a84c12e1abee88b1cc05d9
TEG=${ROOT}/config/tasks/typed_evidence_graph.toml
TEG_SHA256=c1237c32941a2b51c13a8e20741bb6697e541115de6c2c607964724c948bbb92
TARGET=${ROOT}/executions/model-check-899-${SLURM_JOB_ID}-gate-convention
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the gate-convention record.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/record_zero_applicable_gate_convention_cpu.sbatch

# Three read-only output files. Fingerprinted BEFORE.
sha256sum "${W3REC}/ARTIFACTS.json" "${W2GATE}" "${TEG}" > "${STAGE}/readonly_before.sha256"
printf '%s  %s\n' "${W3REC_SHA256}" "${W3REC}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${W2GATE_SHA256}" "${W2GATE}" | sha256sum --check --strict
printf '%s  %s\n' "${TEG_SHA256}" "${TEG}" | sha256sum --check --strict

"${PYTHON}" scripts/record_zero_applicable_gate_convention.py \
  --w3-record "${W3REC}" --w2-gate "${W2GATE}" --typed-evidence-graph "${TEG}" \
  --output "${STAGE}/convention" \
  > "${STAGE}/convention_stdout.json" 2> "${STAGE}/convention_stderr.txt"

sha256sum "${W3REC}/ARTIFACTS.json" "${W2GATE}" "${TEG}" > "${STAGE}/readonly_after.sha256"
cmp -s "${STAGE}/readonly_before.sha256" "${STAGE}/readonly_after.sha256" || { printf 'A read-only artifact changed.\n' >&2; exit 1; }

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum scripts/record_zero_applicable_gate_convention.py \
  slurm/record_zero_applicable_gate_convention_cpu.sbatch \
  "${W2GATE}" "${TEG}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
c = root / "convention"
r = json.loads((c / "zero_applicable_gate_convention.json").read_text(encoding="utf-8"))
if (r["scope"].split(".")[0] != "CAMPAIGN-WIDE"
        or r["the_rule"]["required_third_state"] != "NO_APPLICABLE_CONDITIONS"
        or len(r["the_vacuity_shape_five_instances"]) != 5
        or r["modifies_nothing"] is not True
        or "moved the failure rather than removing it" not in r["the_part_that_matters_most"]["the_lesson"]
        or r["do_not_retrofit"]["instruction"] != "Do NOT open a v3 of the W2 gate to retrofit this rule."):
    raise SystemExit("gate-convention verification differs")
verify_frozen_tree(c)
freeze_tree(root, {"artifact_class": "gate_convention_campaign", "status": "passed"})
verify_frozen_tree(root)
print(json.dumps(r["the_part_that_matters_most"], indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/convention/ARTIFACTS.json"
cat "${TARGET}/readonly_after.sha256"
