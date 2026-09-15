#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-042
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-042-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-042-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v7--a717faf7669b6f26
BASE=${ROOT}/src/masld_bench/evaluators/rna_atac_cobolt_5seed_development.py
WRAPPER=${ROOT}/src/masld_bench/evaluators/rna_atac_cobolt_5seed_v7_development.py
TEST=${ROOT}/tests/unit/test_rna_atac_cobolt_5seed_v7_evaluator.py
TARGET=${ROOT}/executions/rna-atac-cobolt-5seed-v7-evaluator-source-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite evaluator source freeze: %s\n' "${TARGET}" >&2
  exit 1
fi

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python -m py_compile "${BASE}" "${WRAPPER}" "${TEST}"
PYTHONPATH=src python -m unittest \
  tests.unit.test_rna_atac_cobolt_5seed_v7_evaluator \
  tests.unit.test_rna_atac_cobolt_5seed_evaluator \
  tests.unit.test_rna_atac_cobolt_evaluator
PYTHONPATH=src python - "${CANDIDATE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.planner import load_frozen_plan, source_lock

plan = load_frozen_plan(Path(sys.argv[1]))
included = plan["source_lock"].get("included_paths")
if source_lock(Path(plan["package_root"]), included_paths=included) != plan["source_lock"]:
    raise SystemExit("v7 evaluator wrapper altered the frozen campaign source closure")
PY

install -d -m 0750 "${STAGE}/source"
install -m 0640 "${BASE}" "${STAGE}/source/rna_atac_cobolt_5seed_development.py"
install -m 0640 "${WRAPPER}" "${STAGE}/source/rna_atac_cobolt_5seed_v7_development.py"
install -m 0640 "${TEST}" "${STAGE}/source/test_rna_atac_cobolt_5seed_v7_evaluator.py"
PYTHONPATH=src python - "${STAGE}" "${BASE}" "${WRAPPER}" "${TEST}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.hashing import sha256_file

stage = Path(sys.argv[1]).resolve(strict=True)
paths = {
    "base_evaluator": Path(sys.argv[2]).resolve(strict=True),
    "v7_entry_point": Path(sys.argv[3]).resolve(strict=True),
    "unit_test": Path(sys.argv[4]).resolve(strict=True),
}
records = {
    role: {
        "path": path.as_posix(),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    for role, path in paths.items()
}
write_json_exclusive(stage / "source_lock.json", {
    "schema_version": "masld-bench-evaluator-source-lock-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v7",
    "records": records,
    "test_job_id": str(__import__("os").environ["SLURM_JOB_ID"]),
    "outcomes_opened": False,
})
print(freeze_tree(stage, {
    "artifact_class": "evaluator_source_lock",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v7",
    "base_evaluator_sha256": records["base_evaluator"]["sha256"],
    "entry_point_sha256": records["v7_entry_point"]["sha256"],
    "outcomes_opened": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${BASE}" "${WRAPPER}" "${TEST}"
