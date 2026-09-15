#!/usr/bin/env bash

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
ACTIVATION=${ROOT}/executions/gse105127-assay-native-activation-21082444
RRBS_SOURCE=${ROOT}/executions/gse105127-rrbs-source-qc-r2-21066008
TARGET_FASTA=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz
TARGET_FASTA_PLAIN=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
TARGET_GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz

if [[ -z ${GSE105127_CAMPAIGN_ID:-} || ! ${GSE105127_CAMPAIGN_ID} =~ ^[0-9a-f]{16}$ ]]; then
  printf 'GSE105127_CAMPAIGN_ID must be a 16-character lowercase hex lock\n' >&2
  exit 2
fi
CAMPAIGN=${ROOT}/executions/gse105127-production-${GSE105127_CAMPAIGN_ID}
CONTROL=${CAMPAIGN}/control
PLAN=${CAMPAIGN}/plan
REFERENCE=${CAMPAIGN}/reference
RRBS_COLLAPSED=${CAMPAIGN}/rrbs_collapsed
RNA_RAW=${CAMPAIGN}/rna_raw
RSEM_REFERENCE=${CAMPAIGN}/rsem_reference
CROSSWALK=${CAMPAIGN}/cpg_crosswalk
RNA_QUANT=${CAMPAIGN}/rna_quant
RNA_CONSOLIDATED=${CAMPAIGN}/rna_consolidated
ACTIVATED=${CAMPAIGN}/activation

cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${ROOT}/src:${ROOT}"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
"${PYTHON}" - "${CONTROL}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
PY
(
  cd "${ROOT}"
  sha256sum --check --strict "${CONTROL}/source-lock.sha256"
)

replay_complete() {
  local final=$1
  if [[ -f ${final}/COMPLETE ]]; then
    "${PYTHON}" - "${final}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
PY
    printf 'Verified immutable replay: %s\n' "${final}"
    return 0
  fi
  if [[ -e ${final} ]]; then
    printf 'Incomplete immutable target requires manual audit: %s\n' "${final}" >&2
    exit 3
  fi
  return 1
}

begin_attempt() {
  local label=$1
  ATTEMPT=${CAMPAIGN}/attempts/${label}-${SLURM_JOB_ID}${SLURM_ARRAY_TASK_ID:+-${SLURM_ARRAY_TASK_ID}}
  FAILED=${ATTEMPT}.failed
  if [[ -e ${ATTEMPT} || -e ${FAILED} ]]; then
    printf 'Refusing to replace attempt: %s\n' "${ATTEMPT}" >&2
    exit 4
  fi
  mkdir -p "${ATTEMPT}"
  trap 'status=$?; if [[ ${status} -ne 0 && -d ${ATTEMPT} && ! -e ${FAILED} ]]; then mv "${ATTEMPT}" "${FAILED}"; fi; exit ${status}' EXIT
}

publish_stage() {
  local stage=$1
  local final=$2
  "${PYTHON}" - "${stage}" "${final}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import publish_directory_noreplace, verify_frozen_tree

stage, final = map(Path, sys.argv[1:])
verify_frozen_tree(stage)
final.parent.mkdir(parents=True, exist_ok=True)
publish_directory_noreplace(stage, final)
verify_frozen_tree(final)
PY
}

finish_attempt() {
  if [[ -d ${ATTEMPT} ]]; then
    find "${ATTEMPT}" -depth -type d -empty -delete
  fi
  if [[ -e ${ATTEMPT} ]]; then
    printf 'Attempt retained because nonempty work remains: %s\n' "${ATTEMPT}" >&2
    exit 5
  fi
  trap - EXIT
}
