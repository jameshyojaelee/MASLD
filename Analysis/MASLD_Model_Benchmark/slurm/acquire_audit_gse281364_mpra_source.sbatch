#!/usr/bin/env bash
#SBATCH --job-name=masld-gse281364-source
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
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
LOCAL_COUNTS="${PROJECT}/GWAS/finemapping/data/seqfunc_external/hu2025_mpra/raw_counts_v2"
LOCAL_ANALYSIS="${PROJECT}/GWAS/finemapping/results/seqfunc/mpra_benchmark/v2"
BASE="${PROJECT}/.agent/scanpy_env/bin/python"
TARGET="${ROOT}/executions/gse281364-mpra-source-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed GSE281364 source stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE281364 source artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/metadata" "${STAGE}/raw" "${STAGE}/validation"

SERIES_URL='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281364&targ=self&view=full&form=text'
SAMPLE_URL='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8619256&targ=self&view=full&form=text'
ENA_URL='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRX26659333&result=read_run&fields=run_accession,experiment_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout,read_count,base_count&format=tsv'
NATURE_URL='https://www.nature.com/articles/s41588-026-02617-8'
TABLES_URL='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-026-02617-8/MediaObjects/41588_2026_2617_MOESM4_ESM.xlsx'
SOURCE_DATA_URL='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-026-02617-8/MediaObjects/41588_2026_2617_MOESM5_ESM.xlsx'
R1_URL='https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR312/065/SRR31280365/SRR31280365_1.fastq.gz'
R2_URL='https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR312/065/SRR31280365/SRR31280365_2.fastq.gz'

download() {
  local url="$1"
  local output="$2"
  curl --fail --location --silent --show-error \
    --retry 8 --retry-delay 10 --retry-connrefused \
    --continue-at - --output "${output}.partial" "${url}"
  mv "${output}.partial" "${output}"
}

download "${SERIES_URL}" "${STAGE}/metadata/GSE281364.soft.txt"
download "${SAMPLE_URL}" "${STAGE}/metadata/GSM8619256.soft.txt"
download "${ENA_URL}" "${STAGE}/metadata/SRX26659333.ena.tsv"
download "${NATURE_URL}" "${STAGE}/metadata/s41588-026-02617-8.html"
download "${TABLES_URL}" "${STAGE}/metadata/41588_2026_2617_MOESM4_ESM.xlsx"
download "${SOURCE_DATA_URL}" "${STAGE}/metadata/41588_2026_2617_MOESM5_ESM.xlsx"
download "${R1_URL}" "${STAGE}/raw/SRR31280365_1.fastq.gz" &
R1_PID=$!
download "${R2_URL}" "${STAGE}/raw/SRR31280365_2.fastq.gz" &
R2_PID=$!
wait "${R1_PID}"
wait "${R2_PID}"

test "$(stat -c '%s' "${STAGE}/raw/SRR31280365_1.fastq.gz")" -eq 4199316437
test "$(stat -c '%s' "${STAGE}/raw/SRR31280365_2.fastq.gz")" -eq 3462222758
printf '%s  %s\n' \
  e2a206d1a8be823073db430665f514e7 "${STAGE}/raw/SRR31280365_1.fastq.gz" \
  e14106cd95597681759eded6e66c4e02 "${STAGE}/raw/SRR31280365_2.fastq.gz" \
  | md5sum --check --strict
gzip -t "${STAGE}/raw/SRR31280365_1.fastq.gz" &
R1_TEST_PID=$!
gzip -t "${STAGE}/raw/SRR31280365_2.fastq.gz" &
R2_TEST_PID=$!

mapfile -t LOCAL_FILES < <(find "${LOCAL_COUNTS}" -maxdepth 1 -type f -name '*DNA_RNA_counts.txt.gz' | sort)
test "${#LOCAL_FILES[@]}" -eq 16
validation_pids=()
for path in "${LOCAL_FILES[@]}"; do
  (
    expected=$(awk -v target="${path}" '$2 == target {print $1}' "${LOCAL_COUNTS}/sha256sum.txt")
    test "${#expected}" -eq 64
    printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
    gzip -t "${path}"
  ) &
  validation_pids+=("$!")
  if [[ "${#validation_pids[@]}" -eq 4 ]]; then
    for pid in "${validation_pids[@]}"; do wait "${pid}"; done
    validation_pids=()
  fi
done
for pid in "${validation_pids[@]}"; do wait "${pid}"; done
wait "${R1_TEST_PID}"
wait "${R2_TEST_PID}"

"${BASE}" "${ROOT}/scripts/audit_gse281364_mpra_sources.py" \
  --series-soft "${STAGE}/metadata/GSE281364.soft.txt" \
  --sample-soft "${STAGE}/metadata/GSM8619256.soft.txt" \
  --ena-report "${STAGE}/metadata/SRX26659333.ena.tsv" \
  --nature-html "${STAGE}/metadata/s41588-026-02617-8.html" \
  --supplementary-tables "${STAGE}/metadata/41588_2026_2617_MOESM4_ESM.xlsx" \
  --source-data "${STAGE}/metadata/41588_2026_2617_MOESM5_ESM.xlsx" \
  --fastq-r1 "${STAGE}/raw/SRR31280365_1.fastq.gz" \
  --fastq-r2 "${STAGE}/raw/SRR31280365_2.fastq.gz" \
  --local-counts-root "${LOCAL_COUNTS}" \
  --local-analysis-root "${LOCAL_ANALYSIS}" \
  --output "${STAGE}/validation/source_audit.json" \
  > "${STAGE}/validation/source_audit.stdout.txt"

sha256sum "${STAGE}"/metadata/* "${STAGE}"/raw/* > "${STAGE}/source_files.sha256"
sha256sum \
  "${ROOT}/scripts/audit_gse281364_mpra_sources.py" \
  "${ROOT}/slurm/acquire_audit_gse281364_mpra_source.sbatch" \
  "${ROOT}/config/datasets/gse281364.toml" \
  "${LOCAL_COUNTS}/sha256sum.txt" \
  "${LOCAL_ANALYSIS}/gate_verdict.json" \
  > "${STAGE}/source_contracts.sha256"
curl --version > "${STAGE}/validation/curl.version.txt"
"${BASE}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/source_audit.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["experimental_replicates"] != 16
    or receipt["replicates_per_context"] != 4
    or receipt["donor_count"] != 0
    or receipt["sequence_join_state"]
    != "join_unresolved_pending_raw_read_reconstruction"
    or not receipt["sequence_inference_from_interval_prohibited"]
    or receipt["eQTL_or_ieQTL_supervision"]
    or receipt["champion_or_sealed_evaluation"]
):
    raise SystemExit("GSE281364 source audit receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse281364_mpra_source_admission",
        "dataset_id": "gse281364",
        "raw_run": "SRR31280365",
        "experimental_replicates": 16,
        "biological_unit": "experimental_replicate",
        "source_rights": "source_terms_apply_no_raw_redistribution",
        "sequence_join_state": "join_unresolved_pending_raw_read_reconstruction",
        "outcome_role": "exposed_development_MPRA_only",
        "sealed_outcomes_loaded": False,
        "champion_eligible": False,
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 MPRA source admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
