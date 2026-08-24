#!/usr/bin/env bash
#SBATCH --job-name=masld-sequence-reference
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
REFERENCE="${ROOT}/executions/reference-build-21062075"
FASTA="${REFERENCE}/GRCh38.p14.primary_analysis.fa.gz"
FAI="${FASTA}.fai"
GZI="${FASTA}.gzi"
GTF="${REFERENCE}/gencode.v49.primary_analysis.annotation.gtf.gz"
CHROM_SIZES="${REFERENCE}/GRCh38.p14.primary_analysis.chrom.sizes"
TARGET="${ROOT}/executions/sequence-reference-contract-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite sequence reference contract: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
(
  cd "${REFERENCE}"
  sha256sum --check --strict SHA256SUMS
)
printf '%s  %s\n' \
  7b454fee45fd6c019d3fc378428b830a9193e2909bc994ab0eaa59d2ca728ad4 \
  "${FASTA}" | sha256sum --check --strict
printf '%s  %s\n' \
  7c396dc498b1c064852b2e868b7d0fb96327d2dbd5c9fe80f43165bbe9dc0ebb \
  "${FAI}" | sha256sum --check --strict
printf '%s  %s\n' \
  2e419928dd00ac70cf8e3b733b6ae7945daa8427df89cb0daff089cb80ccce8b \
  "${GZI}" | sha256sum --check --strict
printf '%s  %s\n' \
  7da8d8b40ca1194dcbf9af6dc15148b769285fb2cee0ea1f1e2a4aaa92bce65c \
  "${GTF}" | sha256sum --check --strict
printf '%s  %s\n' \
  2c49968ec40d099dfd03d6930c7b224688dadb415c34495f36807ccb3f0253bb \
  "${CHROM_SIZES}" | sha256sum --check --strict
gzip -t "${FASTA}"
gzip -t "${GTF}"

awk '$1 != "chrM" { print }' "${CHROM_SIZES}" \
  > "${STAGE}/GRCh38.p14.sequence_model.chrom.sizes"
test "$(wc -l < "${STAGE}/GRCh38.p14.sequence_model.chrom.sizes")" -eq 24
awk 'BEGIN { for (i=1; i<=22; i++) expected[i]="chr" i; expected[23]="chrX"; expected[24]="chrY" }
     $1 != expected[NR] || $2 !~ /^[0-9]+$/ || $2 <= 0 { failed=1 }
     END { if (NR != 24) failed=1; exit failed }' \
  "${STAGE}/GRCh38.p14.sequence_model.chrom.sizes"

module load SAMtools/1.23.1
module load python/3.11.5-GCCcore-13.2.0
samtools faidx "${FASTA}" chr1:1-20 chrX:1-20 chrY:1-20 \
  > "${STAGE}/validation/reference_boundary_fixture.fa"
python - "${STAGE}/validation/reference_boundary_fixture.fa" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1]).resolve(strict=True)
records = {}
name = None
for line in path.read_text(encoding="utf-8").splitlines():
    if line.startswith(">"):
        name = line[1:]
        records[name] = ""
    elif name is not None:
        records[name] += line.upper()
if set(records) != {"chr1:1-20", "chrX:1-20", "chrY:1-20"}:
    raise SystemExit("reference boundary fixture names differ")
if any(len(value) != 20 or re.fullmatch("[ACGTN]+", value) is None for value in records.values()):
    raise SystemExit("reference boundary fixture sequences differ")
PY

cp "${REFERENCE}/SHA256SUMS" "${STAGE}/upstream-reference.SHA256SUMS"
cp "${REFERENCE}/SOURCE.tsv" "${STAGE}/upstream-reference.SOURCE.tsv"
cat > "${STAGE}/reference.json" <<EOF
{
  "annotation_path": "${GTF}",
  "annotation_sha256": "7da8d8b40ca1194dcbf9af6dc15148b769285fb2cee0ea1f1e2a4aaa92bce65c",
  "bed_coordinate_system": "0_based_half_open",
  "build": "GRCh38.p14",
  "fasta_fai_path": "${FAI}",
  "fasta_fai_sha256": "7c396dc498b1c064852b2e868b7d0fb96327d2dbd5c9fe80f43165bbe9dc0ebb",
  "fasta_gzi_path": "${GZI}",
  "fasta_gzi_sha256": "2e419928dd00ac70cf8e3b733b6ae7945daa8427df89cb0daff089cb80ccce8b",
  "fasta_path": "${FASTA}",
  "fasta_sha256": "7b454fee45fd6c019d3fc378428b830a9193e2909bc994ab0eaa59d2ca728ad4",
  "model_contigs": ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"],
  "model_excluded_contigs": ["chrM"],
  "release": "GENCODE_v49",
  "sequence_job_admission": ["gzip_integrity", "checksum", "build_contig", "ref_allele", "coordinate_convention"]
}
EOF
samtools --version > "${STAGE}/validation/samtools.version.txt"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "sequence_reference_contract",
            "build": "GRCh38.p14",
            "annotation": "GENCODE_v49",
            "fasta_sha256": "7b454fee45fd6c019d3fc378428b830a9193e2909bc994ab0eaa59d2ca728ad4",
            "model_chrom_sizes_sha256": sha256_file(
                root / "GRCh38.p14.sequence_model.chrom.sizes"
            ),
            "model_contigs": 24,
        },
    )
)
PY

mv "${STAGE}" "${TARGET}"
printf 'Sequence reference contract: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
