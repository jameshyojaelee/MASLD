#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-freeze-035
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-freeze-035-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-freeze-035-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
EVALUATOR=${ROOT}/src/masld_bench/evaluators/rna_atac_cobolt_5seed_development.py
TEST=${ROOT}/tests/unit/test_rna_atac_cobolt_5seed_evaluator.py
TARGET=${ROOT}/executions/rna-atac-cobolt-5seed-evaluator-source-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite evaluator source freeze: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/source"
install -m 0640 "${EVALUATOR}" "${STAGE}/source/rna_atac_cobolt_5seed_development.py"
install -m 0640 "${TEST}" "${STAGE}/source/test_rna_atac_cobolt_5seed_evaluator.py"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python - "${STAGE}" "${EVALUATOR}" "${TEST}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.hashing import sha256_file

stage = Path(sys.argv[1]).resolve(strict=True)
evaluator = Path(sys.argv[2]).resolve(strict=True)
test = Path(sys.argv[3]).resolve(strict=True)
records = {
    "evaluator": {
        "path": evaluator.as_posix(),
        "sha256": sha256_file(evaluator),
        "size_bytes": evaluator.stat().st_size,
    },
    "unit_test": {
        "path": test.as_posix(),
        "sha256": sha256_file(test),
        "size_bytes": test.stat().st_size,
    },
}
write_json_exclusive(stage / "source_lock.json", {
    "schema_version": "masld-bench-evaluator-source-lock-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v5",
    "records": records,
    "test_job_ids": ["21069698", "21069701", "21069722"],
    "test_status": "7 tests passed on CPU/nslab after final audit-wording correction",
    "outcomes_opened": False,
})
artifact_hash = freeze_tree(stage, {
    "artifact_class": "evaluator_source_lock",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v5",
    "evaluator_sha256": records["evaluator"]["sha256"],
    "outcomes_opened": False,
})
print(artifact_hash)
PY
mv "${STAGE}" "${TARGET}"
printf 'Evaluator source freeze: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
