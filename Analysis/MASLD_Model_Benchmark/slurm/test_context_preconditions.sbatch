#!/usr/bin/env bash
#SBATCH --job-name=masld-context-preconditions
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python -m unittest -v tests.unit.test_context_preconditions
PYTHONPATH="src:." python - <<'PY'
from pathlib import Path
import tomllib

from masld_bench.context_preconditions import require_precondition_only_campaign

root = Path.cwd()
with (root / "config/campaigns/v1_conditional_model.toml").open("rb") as handle:
    campaign = tomllib.load(handle)
require_precondition_only_campaign(campaign)
if campaign["conditional_model_spec_sha256"] != "UNRESOLVED":
    raise SystemExit("conditional spec unexpectedly resolved")
print("conditional campaign remains blocked and non-executable")
PY
