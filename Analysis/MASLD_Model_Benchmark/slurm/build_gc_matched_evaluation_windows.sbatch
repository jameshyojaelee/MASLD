#!/usr/bin/env bash
#SBATCH --job-name=masld-gc-windows
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
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
DIAGNOSTIC="${ROOT}/executions/gc-background-diagnostic-21063845"
FASTA="${SEQUENCE}/GRCh38.p14.sequence_model.fa"
TARGET="${ROOT}/executions/gc-matched-evaluation-windows-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite GC-matched evaluation windows: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLITS}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 9c70775a4d6542139f2da00aaaf10bd82d6f923cedf7c5012998b9a7569b82bf \
  "${CCRE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 39d75c8e034d4c81b32020b14a9a616b6204215dccde18934450d02440e8f850 \
  "${CCRE}/ENCFF420VPZ.bed.gz" | sha256sum --check --strict
printf '%s  %s\n' 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  "${SEQUENCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ed20dceeb315cffa2bc3de050d846eceefc552451fef741f0dd95b0149d10fce \
  "${FASTA}" | sha256sum --check --strict
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 616cf804a612317067ce30fc4bd15b9c9d72ee5af04bc25e97cfdfced289a694 \
  "${DIAGNOSTIC}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python "${ROOT}/scripts/build_gc_matched_evaluation_windows.py" \
  --split-contract "${SPLITS}" \
  --ccre-bed "${CCRE}/ENCFF420VPZ.bed.gz" \
  --fasta "${FASTA}" \
  --diagnostic-artifacts-sha256 616cf804a612317067ce30fc4bd15b9c9d72ee5af04bc25e97cfdfced289a694 \
  --output "${STAGE}"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_gc_matched_evaluation_windows.py" \
  "${ROOT}/slurm/build_gc_matched_evaluation_windows.sbatch" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
contract = json.loads((root / "contract.json").read_text(encoding="utf-8"))
print(
    freeze_tree(
        root,
        {
            "artifact_class": "genomic_evaluation_windows",
            "positive_windows": contract["positive_windows"],
            "background_windows": contract["background_windows"],
            "background_excludes_all_encode_ccres": True,
            "outcome_independent": True,
            "held_atac_used": False,
            "seed": contract["seed"],
        },
    )
)
PY

mv "${STAGE}" "${TARGET}"
printf 'GC-matched evaluation windows: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
