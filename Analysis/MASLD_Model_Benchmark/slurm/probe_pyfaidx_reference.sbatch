#!/usr/bin/env bash
#SBATCH --job-name=masld-pyfaidx-probe
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
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
OVERLAY="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
REFERENCE="${ROOT}/executions/reference-build-21062075"
CONTRACT="${ROOT}/executions/sequence-reference-contract-21063670"
FASTA="${REFERENCE}/GRCh38.p14.primary_analysis.fa.gz"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d \
  "${CONTRACT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 7b454fee45fd6c019d3fc378428b830a9193e2909bc994ab0eaa59d2ca728ad4 \
  "${FASTA}" | sha256sum --check --strict

module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${FASTA}" <<'PY'
from pathlib import Path
import json
import re
import sys

import pyfaidx

path = Path(sys.argv[1]).resolve(strict=True)
reference = pyfaidx.Fasta(
    str(path),
    as_raw=True,
    sequence_always_upper=True,
    rebuild=False,
    read_ahead=1_000_000,
)
try:
    windows = {
        "chr1": str(reference["chr1"][10_000:12_114]),
        "chrX": str(reference["chrX"][20_000:22_114]),
        "chrY": str(reference["chrY"][30_000:32_114]),
    }
finally:
    reference.close()
if any(len(value) != 2114 for value in windows.values()):
    raise SystemExit("pyfaidx random-access length differs")
if any(re.fullmatch("[ACGTN]+", value) is None for value in windows.values()):
    raise SystemExit("pyfaidx random-access alphabet differs")
complement = str.maketrans("ACGTN", "TGCAN")
sequence = windows["chr1"]
reverse_complement = sequence.translate(complement)[::-1]
if reverse_complement.translate(complement)[::-1] != sequence:
    raise SystemExit("reference reverse-complement fixture failed")
print(
    json.dumps(
        {
            "status": "pass",
            "pyfaidx": pyfaidx.__version__,
            "bgzip_random_access": True,
            "window_length": 2114,
            "reverse_complement_involution": True,
        },
        sort_keys=True,
    )
)
PY
