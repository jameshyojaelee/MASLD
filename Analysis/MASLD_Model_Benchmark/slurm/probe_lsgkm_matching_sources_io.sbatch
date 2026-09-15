#!/usr/bin/env bash
#SBATCH --job-name=model-net-261
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/model-net-261-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite source-connectivity probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

probe() {
  local label="$1"
  local url="$2"
  local output="${STAGE}/${label}.headers.txt"
  local status="${STAGE}/${label}.status.txt"
  set +e
  curl --fail --location --head --proto '=https' --tlsv1.2 \
    --connect-timeout 8 --max-time 15 --silent --show-error \
    "${url}" > "${output}" 2> "${STAGE}/${label}.stderr.txt"
  printf '%s\n' "$?" > "${status}"
  set -e
}

probe bismap_k100_bedgraph \
  https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bedgraph.gz
probe repeatmasker_hg38_dfam2 \
  https://www.repeatmasker.org/genomes/hg38/rmsk4.0.6_dfam2.0/hg38.fa.out.gz

sha256sum "${ROOT}/slurm/probe_lsgkm_matching_sources_io.sbatch" > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
statuses = {
    path.name.removesuffix(".status.txt"): int(path.read_text().strip())
    for path in root.glob("*.status.txt")
}
freeze_tree(root, {
    "artifact_class": "lsgkm_matching_source_connectivity_probe",
    "http_head_exit_codes": statuses,
    "outcomes_read": False,
    "sealed_assets_read": False,
    "biological_training_executed": False,
    "production_predictions_generated": False,
    "status": "diagnostic_only",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Source-connectivity probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
