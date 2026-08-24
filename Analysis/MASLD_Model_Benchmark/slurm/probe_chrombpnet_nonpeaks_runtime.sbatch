#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-nonpeaks-probe
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
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
OVERLAY="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
PREP_OVERLAY="${ROOT}/executions/environments/chrombpnet-native-preprocessing-overlay-21063949"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
BEDTOOLS_RUNTIME="${ROOT}/executions/environments/bedtools-bioconda-2.31.1-21063932"
BEDTOOLS_BIN="${BEDTOOLS_RUNTIME}/runtime/bin"
TARGET="${ROOT}/executions/chrombpnet-nonpeaks-runtime-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${PREP_OVERLAY}/site-packages:${OVERLAY}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite nonpeak runtime probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/mplconfig"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' f1fc72d551449ec0f917d5453183870c8486720f72ff27d963adad5685720f73 \
  "${PREP_OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 04bd59c825874f2d2d7dbf38b741633ccacad4735e95997a7477288e9c9656c4 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 81ed99952bb5caf92cb129e7ea87713468f1eebb143ba0efe6542bfd1a00bc68 \
  "${BEDTOOLS_RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
export SINGULARITYENV_PREPEND_PATH="${BEDTOOLS_BIN}"
"${BEDTOOLS_BIN}/bedtools" --version > "${STAGE}/bedtools.version.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  "${IMAGE}" bedtools --version > "${STAGE}/container-bedtools.version.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env "MPLCONFIGDIR=${STAGE}/mplconfig" \
  --env MPLBACKEND=Agg \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

import pyfaidx

root = Path(sys.argv[1]).resolve(strict=True)
contigs = ("chrTrain", "chrValid", "chrTest")
sequence = "ACGT" * 250
with (root / "fixture.fa").open("x", encoding="utf-8") as handle:
    for contig in contigs:
        handle.write(f">{contig}\n{sequence}\n")
reference = pyfaidx.Fasta(str(root / "fixture.fa"), rebuild=True)
reference.close()
(root / "fixture.chrom.sizes").write_text(
    "".join(f"{contig}\t{len(sequence)}\n" for contig in contigs),
    encoding="utf-8",
)
with (root / "fixture.peaks.bed").open("x", encoding="utf-8") as handle:
    for index, contig in enumerate(contigs):
        start = 200 + 100 * index
        handle.write(
            f"{contig}\t{start}\t{start + 20}\tpeak_{index}\t100\t.\t10\t5\t5\t10\n"
        )
(root / "fixture.blacklist.bed").write_text(
    "chrTrain\t600\t620\n", encoding="utf-8"
)
(root / "fixture.fold.json").write_text(
    json.dumps(
        {"train": ["chrTrain"], "valid": ["chrValid"], "test": ["chrTest"]},
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
PY

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env "MPLCONFIGDIR=${STAGE}/mplconfig" \
  --env MPLBACKEND=Agg \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m chrombpnet.CHROMBPNET prep nonpeaks \
    --genome "${STAGE}/fixture.fa" \
    --output-prefix "${STAGE}/fixture" \
    --peaks "${STAGE}/fixture.peaks.bed" \
    --chrom-sizes "${STAGE}/fixture.chrom.sizes" \
    --chr-fold-path "${STAGE}/fixture.fold.json" \
    --inputlen 20 \
    --stride 10 \
    --neg-to-pos-ratio-train 2 \
    --blacklist-regions "${STAGE}/fixture.blacklist.bed" \
    --seed 20260824 \
    > "${STAGE}/chrombpnet.stdout.txt" \
    2> "${STAGE}/chrombpnet.stderr.txt"

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env "MPLCONFIGDIR=${STAGE}/mplconfig" \
  --env MPLBACKEND=Agg \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
output = root / "fixture_negatives.bed"
rows = [line.rstrip("\n").split("\t") for line in output.open(encoding="utf-8")]
if len(rows) != 5 or any(len(row) != 10 for row in rows):
    raise SystemExit("ChromBPNet nonpeak output shape differs")
fold = {"chrTrain": "train", "chrValid": "valid", "chrTest": "test"}
counts = {value: 0 for value in fold.values()}
for row in rows:
    if row[0] not in fold or int(row[2]) <= int(row[1]) or int(row[9]) != 10:
        raise SystemExit("ChromBPNet nonpeak output geometry differs")
    counts[fold[row[0]]] += 1
if counts != {"train": 2, "valid": 2, "test": 1}:
    raise SystemExit(f"ChromBPNet nonpeak split counts differ: {counts!r}")
result = {
    "schema_version": "masld-bench-chrombpnet-nonpeaks-probe-v1",
    "status": "pass",
    "release": "1.0.1",
    "source_revision": "eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c",
    "bedtools_runtime_artifacts_sha256": "81ed99952bb5caf92cb129e7ea87713468f1eebb143ba0efe6542bfd1a00bc68",
    "preprocessing_overlay_artifacts_sha256": "f1fc72d551449ec0f917d5453183870c8486720f72ff27d963adad5685720f73",
    "input_length": 20,
    "stride": 10,
    "seed": 20260824,
    "output_rows": len(rows),
    "split_counts": counts,
    "blacklist_path_exercised": True,
}
(root / "probe.json").write_text(
    json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(result, sort_keys=True))
PY

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env "MPLCONFIGDIR=${STAGE}/mplconfig" \
  --env MPLBACKEND=Agg \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum "${ROOT}/slurm/probe_chrombpnet_nonpeaks_runtime.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "chrombpnet_nonpeaks_runtime_probe",
            "model_id": "chrombpnet",
            "release": "1.0.1",
            "source_revision": "eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c",
            "bedtools_runtime_artifacts_sha256": "81ed99952bb5caf92cb129e7ea87713468f1eebb143ba0efe6542bfd1a00bc68",
            "preprocessing_overlay_artifacts_sha256": "f1fc72d551449ec0f917d5453183870c8486720f72ff27d963adad5685720f73",
            "seed": 20260824,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet nonpeak runtime probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
