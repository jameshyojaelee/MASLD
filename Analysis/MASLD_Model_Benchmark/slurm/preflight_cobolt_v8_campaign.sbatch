#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-045
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-045-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-045-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/cobolt-v8-campaign-preflight-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite v8 campaign preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/candidates"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config \
  >"${STAGE}/registry_validate.stdout.txt" \
  2>"${STAGE}/registry_validate.stderr.txt"
PYTHONPATH=src python -m masld_bench.cli campaign freeze \
  --campaign config/campaigns/v1_rna_atac_cobolt_5seed_smoke.toml \
  --config-root config \
  --output-root "${STAGE}/candidates" \
  >"${STAGE}/campaign_freeze.stdout.txt" \
  2>"${STAGE}/campaign_freeze.stderr.txt"
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
from collections import Counter
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file
from masld_bench.planner import load_frozen_plan, source_lock

stage = Path(sys.argv[1]).resolve(strict=True)
candidates = sorted(
    stage.joinpath("candidates").glob("v1-rna-atac-cobolt-5seed-smoke-v8--*")
)
if len(candidates) != 1:
    raise SystemExit(f"expected one v8 preflight candidate, observed {len(candidates)}")
candidate = candidates[0]
verify_frozen_tree(candidate)
plan = load_frozen_plan(candidate)
if plan["campaign"]["campaign_id"] != "v1-rna-atac-cobolt-5seed-smoke-v8":
    raise SystemExit("preflight campaign identity differs")
if source_lock(
    Path(plan["package_root"]),
    included_paths=plan["source_lock"].get("included_paths"),
) != plan["source_lock"]:
    raise SystemExit("preflight source lock differs from current tree")
counts = Counter(run["resource_profile"] for run in plan["runs"])
if counts != {"cpu_baseline_smoke": 150, "gpu_rna_atac_smoke": 25}:
    raise SystemExit(f"preflight run census differs: {dict(counts)}")
if any(run.get("task_id") != "rna_conditioned_atac" for run in plan["runs"]):
    raise SystemExit("preflight task binding differs")
if plan["campaign"].get("allow_arrays") is not False:
    raise SystemExit("preflight array policy differs")
write_json_exclusive(stage / "preflight_receipt.json", {
    "schema_version": "masld-bench-cobolt-v8-campaign-preflight-v1",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "preflight_candidate": candidate.relative_to(stage).as_posix(),
    "candidate_artifacts_sha256": sha256_file(candidate / "ARTIFACTS.json"),
    "plan_sha256": plan["plan_sha256"],
    "source_lock_sha256": plan["source_lock"]["source_lock_sha256"],
    "resource_snapshot_sha256": plan["resource_firewall"]["snapshot_sha256"],
    "planned_runs": len(plan["runs"]),
    "cpu_baseline_runs": counts["cpu_baseline_smoke"],
    "gpu_model_runs": counts["gpu_rna_atac_smoke"],
    "arrays": False,
    "runs_submitted": 0,
    "outcomes_opened": False,
    "status": "pass",
})
print(freeze_tree(stage, {
    "artifact_class": "cobolt_v8_campaign_preflight",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "plan_sha256": plan["plan_sha256"],
    "planned_runs": len(plan["runs"]),
    "runs_submitted": 0,
    "outcomes_opened": False,
    "status": "pass",
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
