#!/usr/bin/env bash
#SBATCH --job-name=model-check-895
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
FROZEN=${ROOT}/executions/model-eval-894-21152793-lineage-composition
FROZEN_ARTIFACTS_SHA256=$(sha256sum "${FROZEN}/evaluation/ARTIFACTS.json" | cut -d' ' -f1)
TARGET=${ROOT}/executions/model-check-895-${SLURM_JOB_ID}-closure-diagnostic
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the closure diagnostic.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/diagnose_lineage_closure_subcomposition_cpu.sbatch

# The frozen evaluation is read-only to this step. Fingerprinted before and after.
find "${FROZEN}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_before.txt"

"${PYTHON}" -m unittest tests.unit.test_bulk_lineage_composition_fibrosis   > "${STAGE}/unit_tests.txt" 2>&1

"${PYTHON}" scripts/diagnose_lineage_closure_subcomposition.py \
  --fixture "${FROZEN}/substrate/fixture" \
  --frozen-evaluation "${FROZEN}/evaluation" \
  --output "${STAGE}/diagnostic" \
  --seed 20260826 --replicates 10000 \
  > "${STAGE}/diagnostic_stdout.json" 2> "${STAGE}/diagnostic_stderr.txt"

find "${FROZEN}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_after.txt"
if ! cmp -s "${STAGE}/frozen_before.txt" "${STAGE}/frozen_after.txt"; then
  printf 'The frozen evaluation changed during the diagnostic.\n' >&2; exit 1
fi

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum scripts/diagnose_lineage_closure_subcomposition.py \
  scripts/analyze_bulk_lineage_composition_fibrosis.py \
  slurm/diagnose_lineage_closure_subcomposition_cpu.sbatch \
  "${FROZEN}/evaluation/ARTIFACTS.json" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
d = root / "diagnostic"
r = json.loads((d / "closure_diagnostic_receipt.json").read_text(encoding="utf-8"))
if (r.get("frozen_verdict_altered") is not False
        or r.get("frozen_gate_altered") is not False
        or r.get("is_a_diagnostic_not_a_rescore") is not True
        or not r.get("verdicts")):
    raise SystemExit("closure diagnostic verification differs")
freeze_tree(d, {"artifact_class": "lineage_closure_subcomposition_diagnostic",
                "frozen_verdict_altered": False, "status": "passed"})
verify_frozen_tree(d)
freeze_tree(root, {"artifact_class": "lineage_closure_diagnostic_campaign", "status": "passed"})
verify_frozen_tree(root)
print(json.dumps(r["verdicts"], indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/diagnostic/ARTIFACTS.json"
