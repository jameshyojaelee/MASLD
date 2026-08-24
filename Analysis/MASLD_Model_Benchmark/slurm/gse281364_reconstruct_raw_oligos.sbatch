#!/usr/bin/env bash
#SBATCH --job-name=masld-mpra-reconstruct
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SOURCE="${ROOT}/executions/gse281364-mpra-source-21065555"
SOURCE_SHA=f8e181b00f1297b9b57cee25625d5995c59360f074d2ec3e6caef468fdfd1a20
TOOL="${ROOT}/executions/tools/gse281364-reconstruction-21066228"
TOOL_SHA=52387bdd0a28078608f2342af7154ff962e993c350415f2edc73c3d9c3768307
MAPPING="${PROJECT}/GWAS/finemapping/data/seqfunc_external/hu2025_mpra/raw_counts_v2/GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep1_DNA_RNA_counts.txt.gz"
MAPPING_SHA=7d29c639cc9c59a59b9a74f55421cb2aabc617abe298f657050f688acbd1a57a
R1="${SOURCE}/raw/SRR31280365_1.fastq.gz"
R1_SHA=3e9ecb39576d05594f9719a4777cf9ba16d332b07f75a3661a0189363156c820
R2="${SOURCE}/raw/SRR31280365_2.fastq.gz"
R2_SHA=befba681e707bfd33205d909939058a5b9e264d1f0c5da1d97e41a420ffea3ec
TARGET="${ROOT}/executions/gse281364-raw-reconstruction-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed raw reconstruction: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite raw reconstruction artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/reconstruction" "${STAGE}/validation"
for binding in "${SOURCE_SHA} ${SOURCE}" "${TOOL_SHA} ${TOOL}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
printf '%s  %s\n%s  %s\n%s  %s\n' \
  "${MAPPING_SHA}" "${MAPPING}" \
  "${R1_SHA}" "${R1}" \
  "${R2_SHA}" "${R2}" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" "${TOOL}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("GSE281364 raw reconstruction inputs verified")
PY
module purge

module load GCC/13.2.0
"${TOOL}/bin/gse281364_reconstruct_oligos" \
  --mapping-counts "${MAPPING}" \
  --fastq-r1 "${R1}" \
  --fastq-r2 "${R2}" \
  --output-tsv "${STAGE}/reconstruction/constructs.tsv" \
  --output-stats "${STAGE}/reconstruction/reconstruction_stats.json"
module purge
sha256sum \
  "${ROOT}/slurm/gse281364_reconstruct_raw_oligos.sbatch" \
  "${SOURCE}/ARTIFACTS.json" "${TOOL}/ARTIFACTS.json" "${MAPPING}" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
stats = json.loads((root / "reconstruction/reconstruction_stats.json").read_text())
if (
    stats["status"] != "pass"
    or stats["constructs"] != 10_797
    or stats["first_pass_paired_reads"] != 71_771_807
    or stats["second_pass_paired_reads"] != 71_771_807
    or stats["first_pass_paired_reads"] != stats["second_pass_paired_reads"]
    or stats["observed_count_columns_loaded"]
    or stats["sealed_outcomes_loaded"]
):
    raise SystemExit("raw reconstruction receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse281364_raw_plasmid_oligo_reconstruction",
        "source_artifacts_sha256": "f8e181b00f1297b9b57cee25625d5995c59360f074d2ec3e6caef468fdfd1a20",
        "tool_artifacts_sha256": "52387bdd0a28078608f2342af7154ff962e993c350415f2edc73c3d9c3768307",
        "mapping_source_sha256": "7d29c639cc9c59a59b9a74f55421cb2aabc617abe298f657050f688acbd1a57a",
        "constructs": 10_797,
        "paired_fastq_reads": 71_771_807,
        "sequence_source": "raw_plasmid_paired_fastq_consensus",
        "source_rights": "source_terms_apply_no_raw_or_derived_sequence_redistribution_pending_review",
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "champion_eligible": False,
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 raw reconstruction: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
