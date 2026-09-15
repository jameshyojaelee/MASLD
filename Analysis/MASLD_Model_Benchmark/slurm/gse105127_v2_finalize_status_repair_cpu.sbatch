#!/usr/bin/env bash
#SBATCH --job-name=data-stage-250
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-stage-250-%j.out
#SBATCH --error=executions/logs/data-stage-250-%j.err

source slurm/gse105127_gzip_v2_common.sh
REPAIR_VALIDATION=${ROOT}/executions/data-check-248-21096431
REPAIR_VALIDATION_SHA256=28342689957070508ee8e1a6ca6f2d41cf513a09ce89a5866e42a09de4e825ed
REFERENCE_CONTRACT_SHA256=13903e20b9bb0273a975ee3766b7a2c4f68dcee2cbdbe59608cf4c5010605e63
FINALIZER_REPAIR_SHA256=d02a6604cb80b53d6524f4837b77863b2385698bdd079af2302f5c7ad68c7ada

printf '%s  %s\n' "${REPAIR_VALIDATION_SHA256}" \
  "${REPAIR_VALIDATION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${REFERENCE_CONTRACT_SHA256}" \
  "${ROOT}/scripts/gse105127_exact_reference_contract.py" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FINALIZER_REPAIR_SHA256}" \
  "${ROOT}/scripts/finalize_gse105127_production_activation_status_repair.py" \
  | sha256sum --check --strict
sha256sum --check --strict "${REPAIR_VALIDATION}/source.sha256"
if replay_complete "${TERMINAL}"; then exit 0; fi
begin_attempt finalize-status-repair
if ! replay_complete "${ACTIVATED}"; then
  STAGE=${ATTEMPT}/activation.staging
  "${PYTHON}" scripts/finalize_gse105127_production_activation_status_repair.py \
    --plan-root "${PLAN}" --reference-root "${REFERENCE}" \
    --collapsed-root "${RRBS_COLLAPSED}" --crosswalk-root "${CROSSWALK}" \
    --rna-root "${RNA_CONSOLIDATED}" --output "${STAGE}"
  publish_stage "${STAGE}" "${ACTIVATED}"
fi
STAGE=${ATTEMPT}/terminal.staging
"${PYTHON}" scripts/finalize_gse105127_gzip_preflight_revision_v2.py \
  --control-root "${CONTROL}" --activation-root "${ACTIVATED}" \
  --stage-links-root "${STAGE_LINKS}" \
  --independent-audit-root "${GSE105127_SOURCE_AUDIT_ROOT}" \
  --output "${STAGE}"
publish_stage "${STAGE}" "${TERMINAL}"
finish_attempt
"${PYTHON}" - "${ACTIVATED}" "${TERMINAL}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

activation_root, terminal_root = map(Path, sys.argv[1:])
activation_manifest = verify_frozen_tree(activation_root)
terminal_manifest = verify_frozen_tree(terminal_root)
activation = json.loads(
    (activation_root / "receipt.json").read_text(encoding="utf-8")
)
terminal = json.loads((terminal_root / "receipt.json").read_text(encoding="utf-8"))
if (
    activation_manifest["metadata"].get("artifact_class")
    != "gse105127_production_activation"
    or activation.get("status")
    != "assays_ready_no_fit_evaluator_implementation_pending"
    or activation.get("source_reference_status")
    != "passed_exact_legacy_stream"
    or terminal_manifest["metadata"].get("artifact_class")
    != "gse105127_gzip_preflight_revision_v2_terminal"
    or terminal.get("status") != "assays_ready_no_fit_gzip_contract_repaired"
    or terminal.get("labels_accessed") is not False
    or terminal.get("outcomes_accessed") is not False
    or terminal.get("fit_or_score_performed") is not False
    or terminal.get("clean_or_sealed_champion_eligible_by_this_receipt") is not False
):
    raise SystemExit("repaired terminal GSE105127 contract differs")
PY
sha256sum "${ACTIVATED}/ARTIFACTS.json" "${TERMINAL}/ARTIFACTS.json"
