#!/usr/bin/env bash
#SBATCH --job-name=data-check-237
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-check-237-%j.out
#SBATCH --error=executions/logs/data-check-237-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/data-check-237-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE105127 v2 validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

"${PYTHON}" -m unittest -v tests.unit.test_preflight_gse105127_gzip_inputs_v2 \
  > "${STAGE}/unittest.stdout.txt" 2> "${STAGE}/unittest.stderr.txt"
"${PYTHON}" - <<'PY'
from pathlib import Path

for relative in (
    "scripts/preflight_gse105127_gzip_inputs_v2.py",
    "scripts/prepare_gse105127_gzip_preflight_revision_v2.py",
    "scripts/freeze_gse105127_gzip_stage_link_v2.py",
    "scripts/finalize_gse105127_gzip_preflight_revision_v2.py",
    "scripts/quantify_gse105127_rna_bundle_v2.py",
):
    path = Path(relative)
    compile(path.read_text(encoding="utf-8"), path.as_posix(), "exec")
PY
for wrapper in \
  slurm/gse105127_gzip_v2_common.sh \
  slurm/gse105127_v2_prepare_cpu.sbatch \
  slurm/gse105127_v2_rsem_reference_cpu.sbatch \
  slurm/gse105127_v2_rrbs_collapse_io.sbatch \
  slurm/gse105127_v2_rna_quant_cpu.sbatch \
  slurm/gse105127_v2_crosswalk_io.sbatch \
  slurm/gse105127_v2_consolidate_cpu.sbatch \
  slurm/gse105127_v2_finalize_cpu.sbatch \
  slurm/validate_gse105127_gzip_preflight_v2_cpu.sbatch
do
  bash -n "${wrapper}"
done
"${PYTHON}" - <<'PY' > "${STAGE}/static_contract.stdout.txt"
from pathlib import Path
import tomllib

from scripts.prepare_gse105127_gzip_preflight_revision_v2 import require_runtime_manifest

root = Path.cwd()
require_runtime_manifest(
    root / "executions/environments/transcriptformer-torch2.5.1-21070738"
)
with (root / "config/campaigns/gse105127_gzip_preflight_revision_v2.toml").open("rb") as handle:
    config = tomllib.load(handle)
if (
    config.get("status") != "prepared_not_submitted"
    or config.get("participants") != 19
    or config.get("participant_zone_rows") != 57
    or config.get("strict_gzip", {}).get("member_count") != 1
    or config.get("strict_gzip", {}).get("preflight_runs_in_same_slurm_job_as_consumer") is not True
    or config.get("terminal", {}).get("champion_eligibility") is not False
):
    raise SystemExit("v2 campaign configuration differs")

contracts = {
    "slurm/gse105127_v2_rsem_reference_cpu.sbatch": ("run_preflight rsem_reference", "scripts/build_gse105127_rsem_reference.py"),
    "slurm/gse105127_v2_rrbs_collapse_io.sbatch": ("run_preflight rrbs_collapse", "scripts/collapse_gse105127_rrbs_bundle.py"),
    "slurm/gse105127_v2_rna_quant_cpu.sbatch": ("run_preflight rna_quantification", "scripts/quantify_gse105127_rna_bundle_v2.py"),
    "slurm/gse105127_v2_crosswalk_io.sbatch": ("run_preflight cpg_crosswalk", "scripts/build_gse105127_cpg_crosswalk.py"),
}
for relative, (preflight, consumer) in contracts.items():
    source = (root / relative).read_text(encoding="utf-8")
    if source.index(preflight) >= source.index(consumer) or "freeze_stage_link" not in source:
        raise SystemExit(f"same-job preflight order differs: {relative}")
print("static_campaign_contract=passed")
print("legacy_runtime_checksum_manifest=passed")
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  config/campaigns/gse105127_gzip_preflight_revision_v2.toml \
  scripts/preflight_gse105127_gzip_inputs_v2.py \
  scripts/prepare_gse105127_gzip_preflight_revision_v2.py \
  scripts/freeze_gse105127_gzip_stage_link_v2.py \
  scripts/finalize_gse105127_gzip_preflight_revision_v2.py \
  scripts/quantify_gse105127_rna_bundle_v2.py \
  src/masld_bench/artifacts.py \
  src/masld_bench/hashing.py \
  tests/unit/test_preflight_gse105127_gzip_inputs_v2.py \
  slurm/gse105127_gzip_v2_common.sh \
  slurm/gse105127_v2_prepare_cpu.sbatch \
  slurm/gse105127_v2_rsem_reference_cpu.sbatch \
  slurm/gse105127_v2_rrbs_collapse_io.sbatch \
  slurm/gse105127_v2_rna_quant_cpu.sbatch \
  slurm/gse105127_v2_crosswalk_io.sbatch \
  slurm/gse105127_v2_consolidate_cpu.sbatch \
  slurm/gse105127_v2_finalize_cpu.sbatch \
  slurm/validate_gse105127_gzip_preflight_v2_cpu.sbatch \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, write_json_exclusive

root = Path(sys.argv[1])
tests = (root / "unittest.stderr.txt").read_text(encoding="utf-8")
if "Ran 9 tests" not in tests or "OK" not in tests:
    raise SystemExit("synthetic gzip validation result differs")
write_json_exclusive(
    root / "receipt.json",
    {
        "schema_version": "masld-bench-gse105127-gzip-preflight-v2-validation",
        "status": "passed",
        "synthetic_tests": 9,
        "quantifier_v2_receipt_hash_regression_passed": True,
        "truncated_gzip_rejected": True,
        "trailing_garbage_rejected": True,
        "additional_member_rejected": True,
        "wrong_digest_rejected": True,
        "same_job_wrapper_order_verified": True,
        "molecular_values_opened": False,
        "labels_accessed": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
    },
)
freeze_tree(
    root,
    {
        "artifact_class": "gse105127_gzip_preflight_v2_validation",
        "status": "passed",
        "molecular_values_opened": False,
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
