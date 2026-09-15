#!/usr/bin/env bash
#SBATCH --job-name=model-check-896
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
EVAL=${ROOT}/executions/model-eval-894-21152793-lineage-composition
DIAG=${ROOT}/executions/model-check-895-21153126-closure-diagnostic
GATE_V1=${ROOT}/config/evaluation/bulk_lineage_composition_fibrosis_promotion_gate.json
GATE_V1_SHA256=e9cb9f9d11d1023450b5a50115f0b1b61b326caf99a84c12e1abee88b1cc05d9
GATE_V2=${ROOT}/config/evaluation/bulk_lineage_composition_fibrosis_promotion_gate_v2.json
GATE_V2_SHA256=214da2757089eb3858d2ffed60e34d3ded11efc7e5cab5aa7e946648be1a7674
TASKSPEC=${ROOT}/config/evaluation/bulk_lineage_composition_fibrosis_task.toml
TASKSPEC_SHA256=42f7b71a7f145d10cf2fcddc8b41da015db2f699361a5ad9ce49edc5622ad5b0
TARGET=${ROOT}/executions/model-check-896-${SLURM_JOB_ID}-interpretive-record
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the interpretive record.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/record_lineage_closure_interpretation_cpu.sbatch
printf '%s  %s\n' "${GATE_V1_SHA256}" "${GATE_V1}"   | sha256sum --check --strict
printf '%s  %s\n' "${GATE_V2_SHA256}" "${GATE_V2}"   | sha256sum --check --strict
printf '%s  %s\n' "${TASKSPEC_SHA256}" "${TASKSPEC}" | sha256sum --check --strict

# The frozen evaluation and diagnostic are read-only here. Fingerprinted both ways.
find "${EVAL}" "${DIAG}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_before.txt"

"${PYTHON}" scripts/record_lineage_closure_interpretation.py \
  --evaluation "${EVAL}/evaluation" \
  --diagnostic "${DIAG}/diagnostic" \
  --gate "${GATE_V1}" \
  --taskspec "${TASKSPEC}" \
  --output "${STAGE}/interpretation" \
  > "${STAGE}/record_stdout.json" 2> "${STAGE}/record_stderr.txt"

find "${EVAL}" "${DIAG}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_after.txt"
if ! cmp -s "${STAGE}/frozen_before.txt" "${STAGE}/frozen_after.txt"; then
  printf 'A frozen artifact changed while writing the interpretive record.\n' >&2; exit 1
fi
printf '%s  %s\n' "${GATE_V1_SHA256}" "${GATE_V1}" | sha256sum --check --strict

cp -p "${GATE_V2}" "${STAGE}/forward_only_gate_v2.json"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum scripts/record_lineage_closure_interpretation.py \
  slurm/record_lineage_closure_interpretation_cpu.sbatch \
  "${GATE_V1}" "${GATE_V2}" "${TASKSPEC}" \
  "${EVAL}/evaluation/ARTIFACTS.json" "${DIAG}/diagnostic/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
i = root / "interpretation"
r = json.loads((i / "interpretive_record.json").read_text(encoding="utf-8"))
v2 = json.loads((root / "forward_only_gate_v2.json").read_text(encoding="utf-8"))
if (r.get("frozen_verdict_altered") is not False
        or r.get("frozen_gate_altered") is not False
        or r.get("the_verdict_as_computed", {}).get("verdict") != "PASS"
        or r["the_interpretation"]["qualifying_lineages_after_interpretation"] != ["Hepatocytes"]
        or r["what_the_detection_floor_cost"]["excluded"] != 13
        or r["hepatocyte_depletion_is_a_positive_control"]["this_is_a_statement_about_the_operator_not_a_discovery"] is not True
        or "magnitude_ratio_ignoring_sign" not in json.dumps(r)
        or v2.get("status") != "registered_forward_only_applies_to_future_lanes"
        or v2["does_not_retrofit"]["prior_gate_unchanged"] is not True
        or len(v2["third_instance_of_one_failure_shape"]["instances"]) < 3):
    raise SystemExit("interpretive record verification differs")
verify_frozen_tree(i)
freeze_tree(root, {"artifact_class": "lineage_closure_interpretation_campaign",
                   "frozen_verdict_altered": False, "status": "passed"})
verify_frozen_tree(root)
print(json.dumps(r["the_interpretation"], indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/interpretation/ARTIFACTS.json"
