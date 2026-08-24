#!/usr/bin/env bash
#SBATCH --job-name=masld-sequence-fasta
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE_REFERENCE="${ROOT}/executions/reference-build-21062075"
REFERENCE_CONTRACT="${ROOT}/executions/sequence-reference-contract-21063670"
SOURCE_FASTA="${SOURCE_REFERENCE}/GRCh38.p14.primary_analysis.fa.gz"
SOURCE_FAI="${SOURCE_FASTA}.fai"
MODEL_SIZES="${REFERENCE_CONTRACT}/GRCh38.p14.sequence_model.chrom.sizes"
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
OVERLAY="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
TARGET="${ROOT}/executions/sequence-fasta-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FASTA="${STAGE}/GRCh38.p14.sequence_model.fa"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite sequence FASTA: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 7b454fee45fd6c019d3fc378428b830a9193e2909bc994ab0eaa59d2ca728ad4 \
  "${SOURCE_FASTA}" | sha256sum --check --strict
printf '%s  %s\n' 7c396dc498b1c064852b2e868b7d0fb96327d2dbd5c9fe80f43165bbe9dc0ebb \
  "${SOURCE_FAI}" | sha256sum --check --strict
printf '%s  %s\n' 7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d \
  "${REFERENCE_CONTRACT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 1901777f65ef618569eabb62967c5ace55a63b4bef8795af7a8838fc7042858b \
  "${MODEL_SIZES}" | sha256sum --check --strict
gzip -t "${SOURCE_FASTA}"
gzip -dc "${SOURCE_FASTA}" > "${FASTA}"

module load SAMtools/1.23.1
samtools faidx "${FASTA}"
cut -f1,2 "${FASTA}.fai" > "${STAGE}/validation/materialized.chrom.sizes"
cut -f1,2 "${SOURCE_FAI}" > "${STAGE}/validation/source.chrom.sizes"
cmp \
  "${STAGE}/validation/source.chrom.sizes" \
  "${STAGE}/validation/materialized.chrom.sizes"
awk '$1 != "chrM" { print }' "${STAGE}/validation/materialized.chrom.sizes" \
  | cmp - "${MODEL_SIZES}"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict
module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${FASTA}" \
  > "${STAGE}/validation/pyfaidx-probe.json" <<'PY'
from pathlib import Path
import json
import re
import sys

import pyfaidx

path = Path(sys.argv[1]).resolve(strict=True)
reference = pyfaidx.Fasta(
    str(path),
    as_raw=True,
    sequence_always_upper=True,
    rebuild=False,
    read_ahead=1_000_000,
)
try:
    windows = {
        "chr1": str(reference["chr1"][10_000:12_114]),
        "chrX": str(reference["chrX"][20_000:22_114]),
        "chrY": str(reference["chrY"][30_000:32_114]),
    }
finally:
    reference.close()
if any(len(value) != 2114 for value in windows.values()):
    raise SystemExit("pyfaidx random-access length differs")
if any(re.fullmatch("[ACGTN]+", value) is None for value in windows.values()):
    raise SystemExit("pyfaidx random-access alphabet differs")
complement = str.maketrans("ACGTN", "TGCAN")
sequence = windows["chr1"]
if sequence.translate(complement)[::-1].translate(complement)[::-1] != sequence:
    raise SystemExit("reference reverse-complement fixture failed")
print(
    json.dumps(
        {
            "status": "pass",
            "pyfaidx": pyfaidx.__version__,
            "uncompressed_random_access": True,
            "window_length": 2114,
            "reverse_complement_involution": True,
        },
        sort_keys=True,
    )
)
PY

sha256sum "${FASTA}" "${FASTA}.fai" > "${STAGE}/validation/materialized.sha256"
cat > "${STAGE}/SOURCE.tsv" <<EOF
field	value
build	GRCh38.p14
source_fasta	${SOURCE_FASTA}
source_fasta_sha256	7b454fee45fd6c019d3fc378428b830a9193e2909bc994ab0eaa59d2ca728ad4
materialization	bgzip_decompression_without_sequence_change
model_contigs	chr1-chr22,chrX,chrY
reference_contains_nonmodel_chrM	true
EOF
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
fasta = root / "GRCh38.p14.sequence_model.fa"
print(
    freeze_tree(
        root,
        {
            "artifact_class": "sequence_fasta",
            "build": "GRCh38.p14",
            "source_fasta_sha256": "7b454fee45fd6c019d3fc378428b830a9193e2909bc994ab0eaa59d2ca728ad4",
            "materialized_fasta_sha256": sha256_file(fasta),
            "indexed": True,
            "pyfaidx_random_access_passed": True,
        },
    )
)
PY

mv "${STAGE}" "${TARGET}"
printf 'Sequence FASTA: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
