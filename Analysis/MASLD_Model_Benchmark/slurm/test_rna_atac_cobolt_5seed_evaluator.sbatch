#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-035
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-035-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-035-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v5--162d849f9c6223d5
SOURCE_LOCK=${ROOT}/executions/rna-atac-cobolt-5seed-evaluator-source-21069723
DRIVER=${ROOT}/slurm/evaluate_rna_atac_cobolt_5seed_smoke.sbatch
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
bash -n "${DRIVER}"
PYTHONPATH=src python -m py_compile \
  src/masld_bench/evaluators/rna_atac_cobolt_5seed_development.py \
  tests/unit/test_rna_atac_cobolt_5seed_evaluator.py
PYTHONPATH=src python -m unittest \
  tests.unit.test_rna_atac_cobolt_5seed_evaluator \
  tests.unit.test_rna_atac_cobolt_evaluator
PYTHONPATH=src python - "${CANDIDATE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.planner import load_frozen_plan, source_lock

plan = load_frozen_plan(Path(sys.argv[1]))
included = plan["source_lock"].get("included_paths")
observed = source_lock(Path(plan["package_root"]), included_paths=included)
if observed != plan["source_lock"]:
    raise SystemExit("new five-seed evaluator altered the frozen v5 source closure")
print("five-seed evaluator unit and scoped-source-lock checks passed")
PY
PYTHONPATH=src python - "${SOURCE_LOCK}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
print("five-seed evaluator source freeze verified")
PY
printf '%s  %s\n' \
  6165ece8a32c32e2a636eaae488633d43b43b5fadc981d3f7960181d7c587b88 "${SOURCE_LOCK}/ARTIFACTS.json" \
  43ce6d4d2c08e3191feb6860d73173593586d3febb343b12904f814090fbe8a1 src/masld_bench/evaluators/rna_atac_cobolt_5seed_development.py \
  | sha256sum --check --strict
sha256sum \
  src/masld_bench/evaluators/rna_atac_cobolt_5seed_development.py \
  tests/unit/test_rna_atac_cobolt_5seed_evaluator.py
