#!/usr/bin/env bash
#SBATCH --job-name=model-data-278
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
REFERENCE="${ROOT}/executions/sequence-reference-contract-21063670"
TARGET="${ROOT}/executions/model-data-278-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed matching-track acquisition: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite matching-track acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d \
  "${REFERENCE}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python "${ROOT}/scripts/acquire_lsgkm_matching_tracks.py" \
    --chrom-sizes "${REFERENCE}/GRCh38.p14.sequence_model.chrom.sizes" \
    --output "${STAGE}/tracks" \
    > "${STAGE}/acquisition.stdout.json" \
    2> "${STAGE}/acquisition.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
singularity --version > "${STAGE}/singularity.version.txt"
curl --version > "${STAGE}/curl.version.txt"
sha256sum \
  "${ROOT}/scripts/acquire_lsgkm_matching_tracks.py" \
  "${ROOT}/slurm/acquire_lsgkm_matching_tracks_io.sbatch" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "tracks/receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_source_acquisition_and_structural_audit"
    or not receipt["umap_validation"]["primary_contig_bounds_match_GRCh38p14"]
    or receipt["umap_validation"]["missing_interval_value"] != 0.0
    or receipt["outcomes_read"]
    or receipt["sealed_assets_read"]
    or receipt["biological_training_executed"]
    or receipt["production_predictions_generated"]
):
    raise SystemExit("matching-track acquisition receipt differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_matching_track_acquisition",
    "mappability_source": "Umap_v1.1.0_hg38_k100_multi_read",
    "repeat_source": "RepeatMasker_open_4.0.6_Dfam_2.0_hg38",
    "outcomes_read": False,
    "sealed_assets_read": False,
    "production_training_executed": False,
    "production_predictions_generated": False,
    "status": "passed_source_acquisition",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Matching-track acquisition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
