#!/usr/bin/env bash
#SBATCH --job-name=data-stage-249
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --output=executions/logs/data-stage-249-%j.out
#SBATCH --error=executions/logs/data-stage-249-%j.err

source slurm/gse105127_gzip_v2_common.sh
REPAIR_VALIDATION=${ROOT}/executions/data-check-248-21096431
REPAIR_VALIDATION_SHA256=28342689957070508ee8e1a6ca6f2d41cf513a09ce89a5866e42a09de4e825ed
REFERENCE_CONTRACT_SHA256=13903e20b9bb0273a975ee3766b7a2c4f68dcee2cbdbe59608cf4c5010605e63
CROSSWALK_REPAIR_SHA256=83a73f781eff15e8f3b3849b5496515ea6582b898a8c5828c44e111eb5a2c471

printf '%s  %s\n' "${REPAIR_VALIDATION_SHA256}" \
  "${REPAIR_VALIDATION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${REFERENCE_CONTRACT_SHA256}" \
  "${ROOT}/scripts/gse105127_exact_reference_contract.py" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CROSSWALK_REPAIR_SHA256}" \
  "${ROOT}/scripts/build_gse105127_cpg_crosswalk_status_repair.py" \
  | sha256sum --check --strict
sha256sum --check --strict "${REPAIR_VALIDATION}/source.sha256"
if replay_complete "${STAGE_LINKS}/cpg_crosswalk"; then exit 0; fi
module purge
module load BEDTools/2.31.0-GCC-12.3.0 crossmap/0.7.0
begin_attempt cpg-crosswalk-status-repair
run_preflight cpg_crosswalk cpg_crosswalk \
  --plan-root "${PLAN}" --reference-root "${REFERENCE}" \
  --collapsed-root "${RRBS_COLLAPSED}"
if ! replay_complete "${CROSSWALK}"; then
  STAGE=${ATTEMPT}/output.staging
  "${PYTHON}" scripts/build_gse105127_cpg_crosswalk_status_repair.py \
    --plan-root "${PLAN}" --reference-root "${REFERENCE}" \
    --collapsed-root "${RRBS_COLLAPSED}" --target-fasta "${TARGET_FASTA_PLAIN}" \
    --output "${STAGE}" --work "${ATTEMPT}/work"
  publish_stage "${STAGE}" "${CROSSWALK}"
fi
"${PYTHON}" - "${CROSSWALK}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1])
manifest = verify_frozen_tree(root)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if (
    manifest["metadata"].get("artifact_class")
    != "gse105127_failure_aware_cpg_crosswalk"
    or manifest["metadata"].get("source_reference_status")
    != "passed_exact_legacy_stream"
    or receipt.get("status") != "passed_failure_aware_roundtrip"
    or receipt.get("source_reference_status") != "passed_exact_legacy_stream"
    or receipt.get("all_source_intervals_retained") is not True
    or receipt.get("labels_accessed") is not False
    or receipt.get("fit_or_score_performed") is not False
):
    raise SystemExit("repaired CpG crosswalk contract differs")
PY
freeze_stage_link cpg_crosswalk "${CROSSWALK}" \
  "${STAGE_LINKS}/cpg_crosswalk"
finish_attempt
sha256sum "${CROSSWALK}/ARTIFACTS.json" \
  "${STAGE_LINKS}/cpg_crosswalk/ARTIFACTS.json"
