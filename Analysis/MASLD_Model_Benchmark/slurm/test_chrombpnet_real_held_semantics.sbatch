#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-held-test
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
INPUT="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-inputs-21064244"
FAILED_STAGE="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-seed20260824-21064433.staging"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python - "${INPUT}" "${FAILED_STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree
from scripts.prepare_chrombpnet_train_only_hyperparams import (
    held_rows,
    read_fold,
    read_narrowpeak,
)

source = Path(sys.argv[1]).resolve(strict=True)
attempt = Path(sys.argv[2]).resolve(strict=True)
verify_frozen_tree(source)
fold = read_fold(source / "regions/chrombpnet.fold.json")
pairs = {
    "peaks": (
        source / "regions/chrombpnet.regions.bed",
        attempt / "prepared/native/train_only_filtered.peaks.bed",
    ),
    "nonpeaks": (
        source / "nonpeaks/chrombpnet_negatives.bed",
        attempt / "prepared/native/train_only_filtered.nonpeaks.bed",
    ),
}
for label, (before_path, after_path) in pairs.items():
    before = held_rows(read_narrowpeak(before_path), fold)
    after = held_rows(read_narrowpeak(after_path), fold)
    if before != after:
        raise SystemExit(f"{label} held rows changed semantically")
    print(f"{label}_held_rows={sum(before.values())}")
print("real held-row semantic fixture passed")
PY
