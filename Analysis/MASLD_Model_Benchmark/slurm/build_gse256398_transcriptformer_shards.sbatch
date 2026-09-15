#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-fixture-256
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
SOURCE=${ROOT}/executions/gse256398-source-h5-21069407
CROSSWALK=${ROOT}/executions/gse256398-gencode49-crosswalk-21069439
MEMBERSHIP=${ROOT}/executions/gse256398-transcriptformer-membership-lock-21073862
MEMBERSHIP_SHA256=c3705df69d0bed89a4dcb37947d446f125de0b90f38b21ea615f53e8e73f001e
CHECKPOINT=${ROOT}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens
TARGET=${ROOT}/executions/gse256398-transcriptformer-shards-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

PROMOTION_PENDING=0
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  # The trap used to be cleared one line before the mv, so the post-move
  # verification ran unguarded and a failure left the tree at the
  # authoritative name. A guard switched off immediately before the part that
  # needs it is worse than no guard: a reader sees "trap" and concludes the
  # script is protected.
  if [[ "${PROMOTION_PENDING}" -eq 1 && -e "${TARGET}" && ! -e "${FAILED}" ]]; then
    mv "${TARGET}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE256398 TranscriptFormer fixture.\n' >&2
  exit 1
fi
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH=${ROOT}:${ROOT}/src
cd "${ROOT}"
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_build_gse256398_transcriptformer_shards
"${RUNTIME}/env/bin/python" scripts/build_gse256398_transcriptformer_shards.py \
  --source "${SOURCE}" \
  --crosswalk "${CROSSWALK}" \
  --membership-lock "${MEMBERSHIP}" \
  --membership-lock-sha256 "${MEMBERSHIP_SHA256}" \
  --vocab "${CHECKPOINT}/vocabs/homo_sapiens_gene.h5" \
  --output "${STAGE}"
"${RUNTIME}/env/bin/python" -m pip freeze --all \
  > "${TARGET}.pip-freeze.txt"
mv "${STAGE}" "${TARGET}"
PROMOTION_PENDING=1
PYTHONPATH=${ROOT}/src "${RUNTIME}/env/bin/python" - "${TARGET}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
PY
PROMOTION_PENDING=0
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
