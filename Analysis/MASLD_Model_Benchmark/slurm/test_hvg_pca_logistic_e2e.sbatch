#!/usr/bin/env bash
#SBATCH --job-name=masld-cell-baseline-e2e
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027

package_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
environment_lock="${package_root}/executions/environment-locks/cpu_baseline_smoke--8a7dfd2a15a635d5/environment.lock.json"
output_root="${package_root}/executions/tests/synthetic-cell-baseline-${SLURM_JOB_ID}"

cd "${package_root}"
/gpfs/commons/home/jameslee/.local/bin/micromamba run \
  --prefix /gpfs/commons/home/jameslee/micromamba/envs/immune-atlas \
  env \
  PYTHONNOUSERSITE=1 \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  PYTHONPATH=src \
  python tests/scientific/test_hvg_pca_logistic_e2e.py \
    --environment-lock "${environment_lock}" \
    --output-root "${output_root}"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python -c '
import json
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
root = Path(sys.argv[1]) / "evaluations"
targets = sorted(path for path in root.iterdir() if path.is_dir())
if len(targets) != 5:
    raise SystemExit(f"expected five frozen evaluations, observed {len(targets)}")
verified = []
for target in targets:
    manifest = verify_frozen_tree(target)
    result = json.loads((target / "evaluation.json").read_text(encoding="utf-8"))
    verified.append({
        "path": target.as_posix(),
        "model_id": result["model_id"],
        "artifact_count": len(manifest["artifacts"]),
    })
if len({item["model_id"] for item in verified}) != 5:
    raise SystemExit("synthetic evaluations do not cover five distinct models")
print(json.dumps({"verified_evaluations": verified}, sort_keys=True))
' "${output_root}"
