#!/usr/bin/env bash

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE_CAMPAIGN=${ROOT}/executions/gse105127-production-9675cd48993468d4
PLAN=${SOURCE_CAMPAIGN}/plan
REFERENCE=${SOURCE_CAMPAIGN}/reference
RNA_RAW=${SOURCE_CAMPAIGN}/rna_raw
TARGET_FASTA=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz
TARGET_FASTA_PLAIN=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
TARGET_GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz

: "${GSE105127_REVISION_ID:?GSE105127_REVISION_ID must identify the frozen v2 control}"
: "${GSE105127_SOURCE_AUDIT_ROOT:?GSE105127_SOURCE_AUDIT_ROOT must name the independent audit}"
: "${GSE105127_SOURCE_AUDIT_SHA256:?GSE105127_SOURCE_AUDIT_SHA256 must pin the independent audit}"
: "${GSE105127_V2_VALIDATION_ROOT:?GSE105127_V2_VALIDATION_ROOT must name the frozen v2 validation}"
: "${GSE105127_V2_VALIDATION_SHA256:?GSE105127_V2_VALIDATION_SHA256 must pin the frozen v2 validation}"
if [[ ! ${GSE105127_REVISION_ID} =~ ^[0-9a-f]{16}$ ]]; then
  printf 'GSE105127_REVISION_ID must be a 16-character lowercase hex lock\n' >&2
  exit 2
fi

CAMPAIGN=${ROOT}/executions/gse105127-gzip-preflight-v2-${GSE105127_REVISION_ID}
CONTROL=${CAMPAIGN}/control
PREFLIGHTS=${CAMPAIGN}/preflights
STAGE_LINKS=${CAMPAIGN}/stage_links
RRBS_COLLAPSED=${CAMPAIGN}/rrbs_collapsed
RSEM_REFERENCE=${CAMPAIGN}/rsem_reference
RNA_QUANT=${CAMPAIGN}/rna_quant
CROSSWALK=${CAMPAIGN}/cpg_crosswalk
RNA_CONSOLIDATED=${CAMPAIGN}/rna_consolidated
ACTIVATED=${CAMPAIGN}/activation
TERMINAL=${CAMPAIGN}/terminal

cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${ROOT}/src:${ROOT}"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${GSE105127_SOURCE_AUDIT_SHA256}" \
  "${GSE105127_SOURCE_AUDIT_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${GSE105127_V2_VALIDATION_SHA256}" \
  "${GSE105127_V2_VALIDATION_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
"${PYTHON}" - "${CONTROL}" "${GSE105127_REVISION_ID}" \
  "${GSE105127_SOURCE_AUDIT_ROOT}" "${GSE105127_SOURCE_AUDIT_SHA256}" \
  "${GSE105127_V2_VALIDATION_ROOT}" "${GSE105127_V2_VALIDATION_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

control, revision_id, audit_root, audit_sha, validation_root, validation_sha = sys.argv[1:]
verify_frozen_tree(Path(control))
receipt = json.loads((Path(control) / "control.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "source_audited_revision_locked_before_submission"
    or receipt.get("revision_id") != revision_id
    or receipt.get("independent_source_audit_root") != str(Path(audit_root).resolve(strict=True))
    or receipt.get("independent_source_audit_artifacts_sha256") != audit_sha
    or receipt.get("validation_root") != str(Path(validation_root).resolve(strict=True))
    or receipt.get("validation_artifacts_sha256") != validation_sha
    or receipt.get("labels_accessed") is not False
    or receipt.get("fit_or_score_performed") is not False
):
    raise SystemExit("GSE105127 v2 control contract differs")
PY
sha256sum --check --strict "${CONTROL}/source-lock.sha256"

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
  local job_tag=${SLURM_JOB_ID}${SLURM_ARRAY_TASK_ID:+-${SLURM_ARRAY_TASK_ID}}
  ATTEMPT=${CAMPAIGN}/attempts/${label}-${job_tag}
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

run_preflight() {
  local stage_name=$1
  local logical_name=$2
  shift 2
  local job_tag=${SLURM_JOB_ID}${SLURM_ARRAY_TASK_ID:+-${SLURM_ARRAY_TASK_ID}}
  local staging=${ATTEMPT}/preflight.staging
  PREFLIGHT_FINAL=${PREFLIGHTS}/${logical_name}/${job_tag}
  if [[ -e ${PREFLIGHT_FINAL} ]]; then
    printf 'Refusing to reuse a preflight from this immutable job identity: %s\n' "${PREFLIGHT_FINAL}" >&2
    exit 5
  fi
  "${PYTHON}" scripts/preflight_gse105127_gzip_inputs_v2.py \
    --stage "${stage_name}" --output "${staging}" "$@"
  publish_stage "${staging}" "${PREFLIGHT_FINAL}"
}

freeze_stage_link() {
  local stage_name=$1
  local stage_root=$2
  local final=$3
  shift 3
  local staging=${ATTEMPT}/stage-link.staging
  "${PYTHON}" scripts/freeze_gse105127_gzip_stage_link_v2.py \
    --stage "${stage_name}" --plan-root "${PLAN}" \
    --preflight-root "${PREFLIGHT_FINAL}" --stage-root "${stage_root}" \
    --output "${staging}" "$@"
  publish_stage "${staging}" "${final}"
}

finish_attempt() {
  if [[ -d ${ATTEMPT} ]]; then
    find "${ATTEMPT}" -depth -type d -empty -delete
  fi
  if [[ -e ${ATTEMPT} ]]; then
    printf 'Attempt retained because nonempty work remains: %s\n' "${ATTEMPT}" >&2
    exit 6
  fi
  trap - EXIT
}
