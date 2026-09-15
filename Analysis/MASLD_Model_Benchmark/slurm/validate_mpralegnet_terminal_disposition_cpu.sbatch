#!/usr/bin/env bash
#SBATCH --job-name=model-check-354
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
DISPOSITION="${ROOT}/config/artifacts/models/mpralegnet/terminal_disposition_20260825.json"
RECEIPT="${ROOT}/executions/model-check-230-21092108"
TARGET="${ROOT}/executions/mpralegnet-terminal-disposition-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"; FAILED="${TARGET}.failed"
preserve() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve EXIT
[[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]] && { echo refusing >&2; exit 1; }
install -d -m 0750 "${STAGE}"
sha256sum --check --strict <<EOF
b24475a8393ed2fd2f23acb65813511d1b2504ddaab666bcd961b7c36ef76f46  ${DISPOSITION}
cdf356b81aa19127f2eff45a8d213a4a9d360b9fed5f81edb4c7972d0cb9874d  ${RECEIPT}/ARTIFACTS.json
EOF
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
python -m unittest -v tests.unit.test_mpralegnet_terminal_disposition \
  > "${STAGE}/unit_tests.stdout.txt" 2> "${STAGE}/unit_tests.stderr.txt"
python3 -m pip freeze --all > "${STAGE}/pip_freeze.txt" 2>/dev/null || true
sha256sum "${DISPOSITION}" "${RECEIPT}/ARTIFACTS.json" \
  "${ROOT}/slurm/validate_mpralegnet_terminal_disposition_cpu.sbatch" \
  "${ROOT}/tests/unit/test_mpralegnet_terminal_disposition.py" > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src:${ROOT}" python - "${STAGE}" "${DISPOSITION}" "${RECEIPT}" <<'PY'
import json, sys
from pathlib import Path
from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
stage, disp_path, receipt = (Path(a) for a in sys.argv[1:4])
verify_frozen_tree(receipt)
disp = json.loads(disp_path.read_text())
meta = json.loads((receipt / "ARTIFACTS.json").read_text())["metadata"]
checks = {
    "receipt_status_terminal_blocked": meta["status"] == "terminal_blocked",
    "receipt_provenance_equivalence_failed": meta["provenance_equivalence_pass"] is False,
    "receipt_not_champion_eligible": meta["standalone_champion_eligible"] is False,
    "receipt_no_metrics": meta["metrics_calculated"] is False,
    "receipt_no_outcomes": meta["outcomes_read"] is False,
    "receipt_no_unrestricted_deserializer": meta["unrestricted_deserializer_used"] is False,
    "disposition_matches_receipt_status": disp["terminal_receipt"]["status"] == meta["status"],
    "disposition_matches_terminal_reason": disp["terminal_receipt"]["terminal_reason"] == meta["terminal_reason"],
    "disposition_blocks_scoring": disp["binding_consequences"]["scoring_authorized"] is False,
    "disposition_blocks_ranking": disp["binding_consequences"]["may_enter_a_shortlist_or_ranking"] is False,
    "disposition_blocks_champion": disp["binding_consequences"]["open_champion_eligible"] is False,
    "disposition_is_record_only": disp["record_only"] is True and disp["consumed_by_runtime"] is False,
}
failed = [k for k, v in checks.items() if v is not True]
write_json_exclusive(stage / "validation.json", {
    "schema_version": "masld-bench-mpralegnet-terminal-disposition-validation-v1",
    "checks": checks, "failed_checks": failed,
    "status": "passed" if not failed else "failed",
    "outcomes_read": False, "metrics_calculated": False,
})
freeze_tree(stage, {"artifact_class": "mpralegnet_terminal_disposition_validation",
    "status": "passed" if not failed else "failed",
    "outcomes_read": False, "metrics_calculated": False, "sealed_data_read": False})
print(json.dumps({"status": "passed" if not failed else "failed", "failed": failed}))
if failed: raise SystemExit("validation failed: " + ", ".join(failed))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
