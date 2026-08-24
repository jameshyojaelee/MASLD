#!/usr/bin/env bash
#SBATCH --job-name=masld-gse105127-rrbs-probe
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1

module load python/3.11.5-GCCcore-13.2.0
python - <<'PY'
from __future__ import annotations

import gzip
import io
import json
import urllib.request

url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2819nnn/GSM2819620/suppl/GSM2819620_CV_6610_207.MCSv3.20161129.hs37.cpg.filtered.CG.bed.gz"
request = urllib.request.Request(
    url, headers={"User-Agent": "masld-bench-rrbs-schema-probe/1.0"}
)
with urllib.request.urlopen(request, timeout=120) as response:
    declared = response.headers.get("Content-Length")
    with gzip.GzipFile(fileobj=response) as handle:
        rows = [handle.readline().decode("utf-8", errors="strict").rstrip("\n") for _ in range(20)]
print(json.dumps({"url": url, "content_length": declared, "rows": rows}, indent=2))
PY
