#!/usr/bin/env bash
#SBATCH --job-name=model-check-270
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
INPUT="${ROOT}/executions/model-data-267-21096755.failed/tracks/source/k100.umap.bedgraph.gz"
TARGET="${ROOT}/executions/model-check-270-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Umap census: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
python - "${INPUT}" "${STAGE}/census.json" <<'PY'
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

source = Path(sys.argv[1]).resolve(strict=True)
output = Path(sys.argv[2])
digest = hashlib.sha256()
with source.open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)

counts = Counter()
zero_counts = Counter()
covered = Counter()
minimum = {}
maximum = {}
rows = 0
with gzip.open(source, "rt", encoding="ascii", newline="") as handle:
    for line_number, line in enumerate(handle, start=1):
        if not line.strip() or line.startswith(("#", "browser", "track")):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 4:
            raise SystemExit(f"row width differs at {line_number}")
        contig = fields[0]
        start, end, value = int(fields[1]), int(fields[2]), float(fields[3])
        if start < 0 or end <= start or not math.isfinite(value):
            raise SystemExit(f"row differs at {line_number}")
        counts[contig] += 1
        zero_counts[contig] += value == 0.0
        covered[contig] += end - start
        minimum[contig] = min(value, minimum.get(contig, value))
        maximum[contig] = max(value, maximum.get(contig, value))
        rows += 1
receipt = {
    "schema_version": "masld-bench-lsgkm-umap-source-census-v1",
    "source_path": str(source),
    "source_size_bytes": source.stat().st_size,
    "source_sha256": digest.hexdigest(),
    "rows": rows,
    "rows_by_contig": dict(sorted(counts.items())),
    "explicit_zero_rows_by_contig": dict(sorted(zero_counts.items())),
    "covered_bases_by_contig": dict(sorted(covered.items())),
    "minimum_by_contig": dict(sorted(minimum.items())),
    "maximum_by_contig": dict(sorted(maximum.items())),
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
    "artifact_class": "lsgkm_umap_source_census",
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
printf 'Umap source census: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
