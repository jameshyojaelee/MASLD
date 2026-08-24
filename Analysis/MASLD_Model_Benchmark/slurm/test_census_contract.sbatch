#!/usr/bin/env bash
#SBATCH --job-name=masld-census-contract
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export LC_ALL=C

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

expected_bundle_sha256="${MASLD_VALIDATION_BUNDLE_SHA256:?submit with the explicitly approved validation bundle SHA-256}"
if [[ ! "${expected_bundle_sha256}" =~ ^[0-9a-f]{64}$ ]]; then
  echo "MASLD_VALIDATION_BUNDLE_SHA256 must be a lowercase SHA-256" >&2
  exit 2
fi
bundle_symlink="$(find .gitignore README.md pyproject.toml config src tests slurm/test_census_contract.sbatch -type l -print -quit)"
if [[ -n "${bundle_symlink}" ]]; then
  echo "validation bundle contains a symlink: ${bundle_symlink}" >&2
  exit 2
fi
observed_bundle_sha256="$(
  find .gitignore README.md pyproject.toml config src tests slurm/test_census_contract.sbatch \
    -type f ! -path '*/__pycache__/*' ! -name '*.pyc' -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    | sha256sum \
    | awk '{print $1}'
)"
if [[ "${observed_bundle_sha256}" != "${expected_bundle_sha256}" ]]; then
  echo "validation bundle SHA-256 mismatch: expected ${expected_bundle_sha256}, observed ${observed_bundle_sha256}" >&2
  exit 2
fi
printf 'validation_bundle_sha256=%s\n' "${observed_bundle_sha256}"

module load python/3.11.5-GCCcore-13.2.0

PYTHONPATH=src python -m unittest \
  tests.contract.test_registry_files \
  tests.contract.test_gse296875_phenotype \
  tests.unit.test_artifacts \
  tests.unit.test_checkpoint_preflight \
  tests.unit.test_contracts \
  tests.unit.test_conditional_model \
  tests.unit.test_selection \
  tests.unit.test_slurm \
  tests.unit.test_tournament \
  tests.integration.test_campaign \
  -v
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config
PYTHONPATH=src python -m masld_bench.cli reference validate --config-root config
candidate_root="candidates/census-validation-${SLURM_JOB_ID}"
freeze_json="$(PYTHONPATH=src python -m masld_bench.cli campaign freeze \
  --campaign config/campaigns/v1_admission.toml \
  --config-root config \
  --output-root "${candidate_root}")"
printf '%s\n' "${freeze_json}"
candidate_path="$(python -c 'import json,sys; print(json.load(sys.stdin)["candidate"])' <<<"${freeze_json}")"
PYTHONPATH=src python -m masld_bench.cli campaign verify \
  --candidate "${candidate_path}"
python - "${candidate_path}" <<'PY'
import json
from pathlib import Path
import sys

candidate = Path(sys.argv[1])
plan = json.loads((candidate / "plan.json").read_text(encoding="utf-8"))
review = plan.get("admission_review")
if not isinstance(review, dict) or review.get("review_only") is not True:
    raise SystemExit("admission candidate is not explicitly review-only")
if plan.get("runs") != []:
    raise SystemExit("review-only admission emitted executable runs")
dag = plan.get("dag")
if not isinstance(dag, dict) or any(dag.get(key) != [] for key in ("prerequisites", "run_nodes", "edges")):
    raise SystemExit("review-only admission DAG is not empty")
totals = plan.get("resource_totals")
if not isinstance(totals, dict) or totals.get("jobs") != 0:
    raise SystemExit("review-only admission resource totals are nonzero")
model_ids = sorted(item["model_id"] for item in plan.get("model_dispositions", []))
task_ids = sorted(item["task_id"] for item in plan.get("task_dispositions", []))
if model_ids != review.get("model_ids") or len(model_ids) != review.get("model_count"):
    raise SystemExit("admission review does not cover the complete model census")
if task_ids != review.get("task_ids") or len(task_ids) != review.get("task_count"):
    raise SystemExit("admission review does not cover the complete task census")
variant_records = {
    item["model_id"]: item["variant_capability"]
    for item in plan.get("model_dispositions", [])
    if item.get("variant_capability") is not None
}
if len(variant_records) != 40:
    raise SystemExit("admission review does not bind all variant capability records")
primary_variant_ids = sorted(
    model_id
    for model_id, record in variant_records.items()
    if record.get("primary_eligible") is True
)
if primary_variant_ids != [
    "borzoi_ensemble",
    "context_borzoi",
    "corgi_regular",
    "sequence_only",
    "shuffled_context",
]:
    raise SystemExit("variant primary endpoint roster differs from its frozen capability registry")
if any(
    record.get("primary_eligible") is True
    for record in variant_records.values()
    if record.get("requires_observed_target_context") is True
):
    raise SystemExit("an observed-target-context model entered the sealed primary endpoint")
if (candidate / "jobs").exists() or (candidate / "sbatch_headers.txt").read_text(encoding="utf-8"):
    raise SystemExit("review-only admission emitted job scripts or headers")
print(
    f"review_only_admission models={len(model_ids)} tasks={len(task_ids)} "
    f"variant_models={len(variant_records)} variant_primary={len(primary_variant_ids)} "
    "runs=0 jobs=0"
)
PY
