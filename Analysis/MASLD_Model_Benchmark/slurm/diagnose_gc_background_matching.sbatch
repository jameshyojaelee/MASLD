#!/usr/bin/env bash
#SBATCH --job-name=masld-gc-diagnostic
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SPLITS="${ROOT}/executions/sequence-split-contract-v2-21063827"
CCRE="${ROOT}/executions/encode-v4-ccres-21063708"
SEQUENCE="${ROOT}/executions/sequence-fasta-21063808"
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
OVERLAY="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
FASTA="${SEQUENCE}/GRCh38.p14.sequence_model.fa"
TARGET="${ROOT}/executions/gc-background-diagnostic-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite GC diagnostic: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLITS}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 9c70775a4d6542139f2da00aaaf10bd82d6f923cedf7c5012998b9a7569b82bf \
  "${CCRE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  "${SEQUENCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python scripts/diagnose_gc_background_matching.py \
  --split-contract "${SPLITS}" \
  --ccre-bed "${CCRE}/ENCFF420VPZ.bed.gz" \
  --fasta "${FASTA}" \
  --grid-step 5000 \
  --grid-step 2500 \
  --output "${STAGE}"
sha256sum \
  scripts/diagnose_gc_background_matching.py \
  slurm/diagnose_gc_background_matching.sbatch \
  > "${STAGE}/source.sha256"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/pip-freeze.txt"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
diagnostic = json.loads((root / "diagnostic.json").read_text(encoding="utf-8"))
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gc_background_diagnostic",
            "outcome_independent": diagnostic["outcome_independent"],
            "held_atac_used": diagnostic["held_atac_used"],
            "positive_windows": diagnostic["positive_windows"],
            "scenarios": len(diagnostic["scenarios"]),
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GC diagnostic: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
