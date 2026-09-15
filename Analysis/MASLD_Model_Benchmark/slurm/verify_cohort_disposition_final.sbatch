#!/usr/bin/env bash
#SBATCH --job-name=model-check-903
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a

# Both frozen trees must verify after the config edits.
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import verify_frozen_tree
for r in ("executions/pxd051911-activation-readiness-21109461",
          "executions/gse256398-terminal-disposition-21109462"):
    verify_frozen_tree(Path("${ROOT}") / r)
print("both frozen trees verified")
PY

# Every additive record must cite a live digest for the base registry it leaves alone.
PYTHONPATH=${ROOT}/src python - <<'PY'
import hashlib, json
from pathlib import Path
bad = []
for rec in ("config/campaigns/protein_transport_task_registration_20260825.json",
            "config/campaigns/gse256398_terminal_disposition_20260825.json"):
    d = json.loads(Path(rec).read_text())
    assert d["base_registry_unchanged"] is True, rec
    for entry in d["base_registry_files"]:
        assert entry["modified"] is False, rec
        actual = hashlib.sha256(Path(entry["registry_file"]).read_bytes()).hexdigest()
        if actual != entry["registry_file_sha256"]:
            bad.append(f"{rec}: {entry['registry_file']} {actual} != {entry['registry_file_sha256']}")
    for entry in d["bound_evidence"]:
        p = Path(entry["path"])
        if not p.exists():
            bad.append(f"{rec}: missing bound evidence {p}")
            continue
        actual = hashlib.sha256(p.read_bytes()).hexdigest()
        if actual != entry["sha256"]:
            bad.append(f"{rec}: {p} {actual} != {entry['sha256']}")
if bad:
    raise SystemExit("STALE BINDINGS:\n" + "\n".join(bad))
print("all additive-record bindings rederive")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_build_pxd051911_activation_readiness \
  tests.unit.test_build_gse256398_terminal_disposition \
  tests.unit.test_check_contract_suite_baseline \
  tests.contract.test_registry_files \
  tests.contract.test_cross_cohort_plan \
  > "${ROOT}/executions/cohort-disposition-final-${SLURM_JOB_ID}.unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${ROOT}/executions/cohort-disposition-final-${SLURM_JOB_ID}.unittest.txt" \
  --require-module test_build_pxd051911_activation_readiness \
  --require-module test_build_gse256398_terminal_disposition \
  --require-module test_registry_files
echo "FINAL VERIFICATION PASSED"
