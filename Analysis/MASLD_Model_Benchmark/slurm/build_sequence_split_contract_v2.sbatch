#!/usr/bin/env bash
#SBATCH --job-name=masld-sequence-splits-v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CCRE="${ROOT}/executions/encode-v4-ccres-21063708"
REFERENCE="${ROOT}/executions/sequence-reference-contract-21063670"
SEQUENCE="${ROOT}/executions/sequence-fasta-21063808"
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
OVERLAY="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
FASTA="${SEQUENCE}/GRCh38.p14.sequence_model.fa"
TARGET="${ROOT}/executions/sequence-split-contract-v2-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite sequence split contract: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' 9c70775a4d6542139f2da00aaaf10bd82d6f923cedf7c5012998b9a7569b82bf \
  "${CCRE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 39d75c8e034d4c81b32020b14a9a616b6204215dccde18934450d02440e8f850 \
  "${CCRE}/ENCFF420VPZ.bed.gz" | sha256sum --check --strict
printf '%s  %s\n' 7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d \
  "${REFERENCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 1901777f65ef618569eabb62967c5ace55a63b4bef8795af7a8838fc7042858b \
  "${REFERENCE}/GRCh38.p14.sequence_model.chrom.sizes" | sha256sum --check --strict
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

module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python scripts/build_sequence_split_contract.py \
  --ccre-acquisition "${CCRE}" \
  --ccre-artifacts-sha256 9c70775a4d6542139f2da00aaaf10bd82d6f923cedf7c5012998b9a7569b82bf \
  --ccre-bed-sha256 39d75c8e034d4c81b32020b14a9a616b6204215dccde18934450d02440e8f850 \
  --reference-contract "${REFERENCE}" \
  --reference-artifacts-sha256 7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d \
  --chrom-sizes-sha256 1901777f65ef618569eabb62967c5ace55a63b4bef8795af7a8838fc7042858b \
  --sequence-fasta-contract "${SEQUENCE}" \
  --sequence-fasta-artifacts-sha256 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  --fasta-sha256 ed20dceeb315cffa2bc3de050d846eceefc552451fef741f0dd95b0149d10fce \
  --output "${STAGE}"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum \
  scripts/build_sequence_split_contract.py \
  slurm/build_sequence_split_contract_v2.sbatch \
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
            "artifact_class": "sequence_split_contract",
            "schema_version": contract["schema_version"],
            "dataset_id": "gse296875",
            "genomic_folds": contract["genomic_folds"],
            "crossed_outer_splits": contract["crossed_outer_splits"],
            "ccre_windows": contract["ccre_windows"],
            "held_atac_used_for_window_selection": False,
            "sequence_ambiguity_policy": contract["sequence_ambiguity_policy"],
            "selection_seed": contract["selection_seed"],
        },
    )
)
PY

mv "${STAGE}" "${TARGET}"
printf 'Sequence split contract: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
