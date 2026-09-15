#!/usr/bin/env bash
#SBATCH --job-name=model-check-276
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
INPUT="${ROOT}/executions/model-data-275-21096917.failed/tracks/source/hg38.fa.out.gz"
TARGET="${ROOT}/executions/model-check-276-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite RepeatMasker census: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
python - "${INPUT}" "${STAGE}/census.json" <<'PY'
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import sys

source = Path(sys.argv[1]).resolve(strict=True)
output = Path(sys.argv[2])
digest = hashlib.sha256()
with source.open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)

widths = Counter()
samples = defaultdict(list)
primary_counts = Counter()
rows = 0
with gzip.open(source, "rt", encoding="ascii", newline="") as handle:
    for line_number, line in enumerate(handle, start=1):
        fields = line.split()
        if not fields or not fields[0].isdigit():
            continue
        width = len(fields)
        widths[width] += 1
        if len(samples[width]) < 5:
            samples[width].append({"line_number": line_number, "line": line.rstrip("\n")})
        if width >= 7 and fields[4].startswith("chr"):
            primary_counts[fields[4]] += 1
        rows += 1
receipt = {
    "schema_version": "masld-bench-lsgkm-repeatmasker-source-census-v1",
    "source_path": str(source),
    "source_size_bytes": source.stat().st_size,
    "source_sha256": digest.hexdigest(),
    "rows": rows,
    "row_width_counts": {str(key): value for key, value in sorted(widths.items())},
    "row_width_samples": {str(key): value for key, value in sorted(samples.items())},
    "chr_prefixed_rows": dict(sorted(primary_counts.items())),
    "outcomes_read": False,
    "sealed_assets_read": False,
}
output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {
    "artifact_class": "lsgkm_repeatmasker_source_census",
    "outcomes_read": False,
    "sealed_assets_read": False,
    "production_training_executed": False,
    "production_predictions_generated": False,
    "status": "diagnostic_only",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'RepeatMasker source census: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
