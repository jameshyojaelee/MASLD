#!/usr/bin/env bash
#SBATCH --job-name=masld-macs2-probe
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/macs2-runtime-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite MACS2 probe: %s\n' "${TARGET}" >&2
  exit 1
fi
mkdir -m 750 "${STAGE}"
module load MACS2/2.2.9.1-foss-2023b
macs2 --version > "${STAGE}/macs2.version.txt"
python - "${STAGE}/fixture.tagAlign.gz" <<'PY'
from pathlib import Path
import gzip
import sys

path = Path(sys.argv[1])
with gzip.open(path, "xt", encoding="utf-8") as handle:
    for index in range(200):
        position = 4900 + index % 25
        strand = "+" if index % 2 == 0 else "-"
        handle.write(f"chr1\t{position}\t{position + 1}\tread_{index}\t1000\t{strand}\n")
PY
macs2 callpeak \
  --treatment "${STAGE}/fixture.tagAlign.gz" \
  --format BED \
  --gsize 10000 \
  --name fixture \
  --outdir "${STAGE}" \
  --nomodel \
  --shift -75 \
  --extsize 150 \
  --keep-dup all \
  --qvalue 0.05 \
  --call-summits \
  > "${STAGE}/macs2.stdout.txt" \
  2> "${STAGE}/macs2.stderr.txt"
if [[ ! -s "${STAGE}/fixture_peaks.narrowPeak" ]]; then
  printf 'MACS2 fixture did not emit a nonempty narrowPeak file\n' >&2
  exit 1
fi
python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum "${ROOT}/slurm/probe_macs2_runtime.sbatch" > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
rows = sum(1 for _line in (root / "fixture_peaks.narrowPeak").open())
print(
    freeze_tree(
        root,
        {
            "artifact_class": "macs2_runtime_probe",
            "version": "2.2.9.1",
            "fixture_narrowpeaks": rows,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'MACS2 runtime probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
