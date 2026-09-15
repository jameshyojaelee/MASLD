#!/usr/bin/env bash
#SBATCH --job-name=data-check-225
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${INDEPENDENT_VALIDATION_ROOT:?INDEPENDENT_VALIDATION_ROOT must name the frozen validation}"
: "${INDEPENDENT_VALIDATION_SHA256:?INDEPENDENT_VALIDATION_SHA256 must pin that ARTIFACTS.json}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
AUDITOR=scripts/audit_gse105127_source_dag_unexposed.py
AUDITOR_SHA256=5b8698983bbbc2e8f03b363d93610848a16f624bd6c2b54550c303f8fe8b8347
CAMPAIGN=${ROOT}/executions/gse105127-production-9675cd48993468d4
TARGET_FASTA=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz
TARGET_GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
TARGET=${ROOT}/executions/data-check-225-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE105127 independent source audit.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${INDEPENDENT_VALIDATION_SHA256}" "${INDEPENDENT_VALIDATION_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${AUDITOR_SHA256}" "${AUDITOR}" | sha256sum --check --strict
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
bash -n slurm/audit_gse105127_source_dag_unexposed_io.sbatch
{
  printf 'JobID\tJobName\tState\tExitCode\tSubmitLine\n'
  sacct -X \
    -j 21083518,21083520,21083522,21087896,21087897,21087898,21088014,21088066,21088068 \
    -n -P -o JobID,JobName,State,ExitCode,SubmitLine%1000 | tr '|' '\t'
} > "${STAGE}/slurm_jobs.tsv"
"${PYTHON}" "${AUDITOR}" \
  --project-root "${ROOT}" \
  --campaign-root "${CAMPAIGN}" \
  --legacy-diagnostic-root "${ROOT}/executions/data-check-211-21084942" \
  --alphabet-diagnostic-root "${ROOT}/executions/data-check-216-21086818" \
  --source-cache-root "${ROOT}/executions/gse105127-reference-source-cache-21087810" \
  --target-fasta "${TARGET_FASTA}" \
  --target-gtf "${TARGET_GTF}" \
  --submission-record "${STAGE}/slurm_jobs.tsv" \
  --output "${STAGE}/audit" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  "${AUDITOR}" \
  slurm/audit_gse105127_source_dag_unexposed_io.sbatch \
  "${INDEPENDENT_VALIDATION_ROOT}/ARTIFACTS.json" \
  "${CAMPAIGN}/control/ARTIFACTS.json" \
  "${CAMPAIGN}/plan/ARTIFACTS.json" \
  "${CAMPAIGN}/reference/ARTIFACTS.json" \
  "${CAMPAIGN}/rsem_reference/ARTIFACTS.json" \
  "${STAGE}/slurm_jobs.tsv" \
  "${STAGE}/audit/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${INDEPENDENT_VALIDATION_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
validation_sha256 = sys.argv[2]
manifest = verify_frozen_tree(root / "audit")
receipt = json.loads((root / "audit/independent_source_dag_audit.json").read_text())
if (
    manifest["metadata"].get("compressed_legacy_source_usable_by_sequence_jobs") is not False
    or manifest["metadata"].get("decompressed_legacy_reference_eligible_for_lift_mechanics") is not True
    or manifest["metadata"].get("strict_every_sequence_job_begins_with_gzip_t_contract_passed") is not False
    or manifest["metadata"].get("current_dag_promotion_eligible") is not False
    or receipt.get("rna_or_methylation_values_parsed") is not False
    or receipt.get("modeling_or_scoring_authorized_by_this_receipt") is not False
    or receipt.get("clean_or_sealed_champion_eligible_by_this_receipt") is not False
):
    raise SystemExit("GSE105127 independent source-only firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse105127_independent_source_dag_audit_execution",
        "validation_artifacts_sha256": validation_sha256,
        "compressed_legacy_source_usable_by_sequence_jobs": False,
        "decompressed_legacy_reference_eligible_for_lift_mechanics": True,
        "strict_every_sequence_job_begins_with_gzip_t_contract_passed": False,
        "current_dag_promotion_eligible": False,
        "molecular_values_opened": False,
        "modeling_or_scoring_authorized": False,
        "clean_or_sealed_champion_eligible": False,
        "status": "passed_audit_source_mechanics_verified_promotion_blocked",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/audit/ARTIFACTS.json"
