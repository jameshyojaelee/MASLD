#!/usr/bin/env bash
#SBATCH --job-name=source-xwalk-256398
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/gse256398-source-h5-21069407
SOURCE_SHA=f54ac799b22cb02912b82432929255e6400d871cc0fd99dbf9b4756f8f977a20
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
GTF_SHA=73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9
TARGET=${ROOT}/executions/gse256398-gencode49-crosswalk-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE256398 crosswalk: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a
PYTHONPATH=${ROOT}:${ROOT}/src:${PYTHONPATH:-} python -m unittest -v \
  tests.unit.test_build_gse256398_gencode49_crosswalk
PYTHONPATH=${ROOT}/src python - "${SOURCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("GSE256398 source frozen-tree verification passed")
PY
gzip -t "${GTF}"
printf '%s  %s\n' "${GTF_SHA}" "${GTF}" | sha256sum --check --strict
PYTHONPATH=${ROOT}:${ROOT}/src:${PYTHONPATH:-} python \
  scripts/build_gse256398_gencode49_crosswalk.py \
  --source "${SOURCE}" \
  --source-artifacts-sha256 "${SOURCE_SHA}" \
  --gencode49-gtf "${GTF}" \
  --gencode49-gtf-sha256 "${GTF_SHA}" \
  --output "${STAGE}"
sha256sum \
  config/datasets/gse256398.toml \
  scripts/build_gse256398_gencode49_crosswalk.py \
  tests/unit/test_build_gse256398_gencode49_crosswalk.py \
  slurm/build_gse256398_gencode49_crosswalk.sbatch \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "summary.json").read_text(encoding="utf-8"))
if (
    summary["status"] != "pass_with_unmapped_features_masked"
    or summary["source_features"] != 36_601
    or summary["model_training_activated"]
    or summary["sealed_or_final_outcomes_read"]
    or not summary["reference_mapping_is_outcome_independent"]
):
    raise SystemExit("GSE256398 crosswalk receipt differs")
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse256398_gencode49_crosswalk",
            "dataset_id": "gse256398",
            "source_features": summary["source_features"],
            "exact_mapped_features": summary["mapping_states"].get(
                "stable_id_exact_gencode_v49", 0
            ),
            "unmapped_features_masked": True,
            "outcome_independent": True,
            "model_training_activated": False,
            "sealed": False,
        },
    )
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE256398 GENCODE v49 crosswalk: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
