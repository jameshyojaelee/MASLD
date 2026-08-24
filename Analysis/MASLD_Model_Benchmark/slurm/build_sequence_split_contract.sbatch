#!/usr/bin/env bash
#SBATCH --job-name=masld-sequence-splits
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
CCRE_ARTIFACTS_SHA256=9c70775a4d6542139f2da00aaaf10bd82d6f923cedf7c5012998b9a7569b82bf
CCRE_BED_SHA256=39d75c8e034d4c81b32020b14a9a616b6204215dccde18934450d02440e8f850
REFERENCE_ARTIFACTS_SHA256=7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d
CHROM_SIZES_SHA256=1901777f65ef618569eabb62967c5ace55a63b4bef8795af7a8838fc7042858b
TARGET="${ROOT}/executions/sequence-split-contract-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite sequence split contract: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' "${CCRE_ARTIFACTS_SHA256}" "${CCRE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CCRE_BED_SHA256}" "${CCRE}/ENCFF420VPZ.bed.gz" \
  | sha256sum --check --strict
printf '%s  %s\n' "${REFERENCE_ARTIFACTS_SHA256}" "${REFERENCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CHROM_SIZES_SHA256}" \
  "${REFERENCE}/GRCh38.p14.sequence_model.chrom.sizes" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python scripts/build_sequence_split_contract.py \
  --ccre-acquisition "${CCRE}" \
  --ccre-artifacts-sha256 "${CCRE_ARTIFACTS_SHA256}" \
  --ccre-bed-sha256 "${CCRE_BED_SHA256}" \
  --reference-contract "${REFERENCE}" \
  --reference-artifacts-sha256 "${REFERENCE_ARTIFACTS_SHA256}" \
  --chrom-sizes-sha256 "${CHROM_SIZES_SHA256}" \
  --output "${STAGE}"
sha256sum \
  scripts/build_sequence_split_contract.py \
  slurm/build_sequence_split_contract.sbatch \
  > "${STAGE}/source.sha256"
python -m pip freeze > "${STAGE}/pip-freeze.txt"
PYTHONPATH="src" python - "${STAGE}" <<'PY'
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
            "dataset_id": "gse296875",
            "genomic_folds": contract["genomic_folds"],
            "crossed_outer_splits": contract["crossed_outer_splits"],
            "ccre_windows": contract["ccre_windows"],
            "held_atac_used_for_window_selection": False,
            "selection_seed": contract["selection_seed"],
        },
    )
)
PY

mv "${STAGE}" "${TARGET}"
printf 'Sequence split contract: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
