#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-046
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-046-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-046-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/cobolt-v8-registry-only-validation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite v8 registry-only validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config \
  >"${STAGE}/registry_validate.stdout.txt" \
  2>"${STAGE}/registry_validate.stderr.txt"
PYTHONPATH=src python -m unittest \
  tests.unit.test_rna_atac_cobolt \
  tests.unit.test_rna_atac_cobolt_evaluator \
  tests.unit.test_rna_atac_cobolt_5seed_evaluator \
  tests.unit.test_rna_atac_cobolt_5seed_v8_evaluator \
  tests.integration.test_campaign.CampaignIntegrationTests.test_scoped_cobolt_lock_ignores_queue_but_detects_bound_source \
  >"${STAGE}/unittest.stdout.txt" \
  2>"${STAGE}/unittest.stderr.txt"
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
import tomllib
import sys

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.hashing import sha256_file

stage = Path(sys.argv[1]).resolve(strict=True)
root = Path.cwd()
campaign_path = root / "config/campaigns/v1_rna_atac_cobolt_5seed_smoke.toml"
campaign = tomllib.loads(campaign_path.read_text(encoding="utf-8"))
if campaign.get("campaign_id") != "v1-rna-atac-cobolt-5seed-smoke-v8":
    raise SystemExit("Cobolt campaign identity differs")
write_json_exclusive(stage / "validation_receipt.json", {
    "schema_version": "masld-bench-cobolt-v8-registry-only-validation-v1",
    "campaign_id": campaign["campaign_id"],
    "campaign_config_sha256": sha256_file(campaign_path),
    "targeted_test_modules": [
        "tests.unit.test_rna_atac_cobolt",
        "tests.unit.test_rna_atac_cobolt_evaluator",
        "tests.unit.test_rna_atac_cobolt_5seed_evaluator",
        "tests.unit.test_rna_atac_cobolt_5seed_v8_evaluator",
        "tests.integration.test_campaign.CampaignIntegrationTests.test_scoped_cobolt_lock_ignores_queue_but_detects_bound_source",
    ],
    "campaign_planned_or_frozen": False,
    "runs_submitted": 0,
    "outcomes_opened": False,
    "status": "pass",
})
print(freeze_tree(stage, {
    "artifact_class": "cobolt_v8_registry_only_validation",
    "campaign_id": campaign["campaign_id"],
    "campaign_planned_or_frozen": False,
    "runs_submitted": 0,
    "outcomes_opened": False,
    "status": "pass",
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
