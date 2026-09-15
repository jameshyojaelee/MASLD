#!/usr/bin/env bash
#SBATCH --job-name=model-check-897
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
EVAL=${ROOT}/executions/model-eval-894-21152793-lineage-composition/evaluation
INTERP=${ROOT}/executions/model-check-896-21153799-interpretive-record/interpretation
GATE=${ROOT}/config/evaluation/bulk_lineage_composition_fibrosis_promotion_gate.json
GATE_SHA256=e9cb9f9d11d1023450b5a50115f0b1b61b326caf99a84c12e1abee88b1cc05d9
TARGET=${ROOT}/executions/model-check-897-${SLURM_JOB_ID}-prereg-evidence
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the evidence record.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/record_w2_preregistration_evidence_cpu.sbatch
printf '%s  %s\n' "${GATE_SHA256}" "${GATE}" | sha256sum --check --strict

find "${EVAL}" "${INTERP}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_before.txt"

"${PYTHON}" scripts/record_w2_preregistration_evidence.py \
  --evaluation "${EVAL}" \
  --gate "${GATE}" \
  --interpretive-record "${INTERP}" \
  --job-id 21152793 \
  --output "${STAGE}/evidence" \
  > "${STAGE}/evidence_stdout.json" 2> "${STAGE}/evidence_stderr.txt"

find "${EVAL}" "${INTERP}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_after.txt"
cmp -s "${STAGE}/frozen_before.txt" "${STAGE}/frozen_after.txt" || { printf 'A frozen artifact changed.\n' >&2; exit 1; }
printf '%s  %s\n' "${GATE_SHA256}" "${GATE}" | sha256sum --check --strict

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum scripts/record_w2_preregistration_evidence.py \
  slurm/record_w2_preregistration_evidence_cpu.sbatch "${GATE}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
e = root / "evidence"
r = json.loads((e / "preregistration_evidence.json").read_text(encoding="utf-8"))
if (r["strong_evidence_no_clock_involved"]["identical"] is not True
        or r.get("record_is_additive_overlay") is not True
        or "caveat" not in r["corroborating_evidence_timestamp_based_and_weaker"]):
    raise SystemExit("pre-registration evidence verification differs")
verify_frozen_tree(e)
freeze_tree(root, {"artifact_class": "w2_preregistration_evidence_campaign", "status": "passed"})
verify_frozen_tree(root)
print(json.dumps(r["strong_evidence_no_clock_involved"], indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/evidence/ARTIFACTS.json"
