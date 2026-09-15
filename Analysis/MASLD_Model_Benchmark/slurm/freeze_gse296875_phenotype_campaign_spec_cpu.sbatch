#!/usr/bin/env bash
#SBATCH --job-name=model-check-001
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/gse296875-phenotype-campaign-spec-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite campaign spec artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHON=python3
mkdir -p "${STAGE}/validation"

# Every frozen input is checked before anything is read from it.
"${PYTHON}" - <<'PY' > "${STAGE}/validation/bound_fixtures.sha256"
from hashlib import sha256
from pathlib import Path
import sys
sys.path.insert(0, "scripts")
import importlib.util
spec = importlib.util.spec_from_file_location(
    "freeze_campaign_spec", "scripts/freeze_gse296875_phenotype_campaign_spec.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
for name, (relative, expected) in sorted(module.BOUND_FIXTURES.items()):
    print(f"{expected}  {relative}/ARTIFACTS.json")
PY
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

PYTHONPATH=src "${PYTHON}" -m unittest -v \
  tests.unit.test_gse296875_phenotype_campaign_spec \
  tests.unit.test_gse296875_phenotype_scoring \
  tests.contract.test_gse296875_phenotype \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH=src "${PYTHON}" scripts/freeze_gse296875_phenotype_campaign_spec.py \
  --root "${ROOT}" \
  --output "${STAGE}/spec" \
  > "${STAGE}/validation/spec_stdout.json"

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true

sha256sum \
  scripts/freeze_gse296875_phenotype_campaign_spec.py \
  src/masld_bench/gse296875_phenotype_scoring.py \
  tests/unit/test_gse296875_phenotype_campaign_spec.py \
  tests/unit/test_gse296875_phenotype_scoring.py \
  tests/contract/test_gse296875_phenotype.py \
  config/evaluation/gse296875_histopathology.toml \
  config/datasets/gse296875.toml \
  slurm/freeze_gse296875_phenotype_campaign_spec_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH=src "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
spec = json.loads((root / "spec" / "campaign_spec.json").read_text())
if spec["champion_eligible"] or spec["role"] != "secondary_development_only":
    raise SystemExit("campaign spec role differs")
if spec["multiple_testing"]["family_size"] != 12:
    raise SystemExit("confirmatory family size differs")
if spec["arms"]["race_is_not_fitted"] is not True:
    raise SystemExit("race policy differs")
if spec["lineages"]["primary"] != [
    "cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell"
]:
    raise SystemExit("frozen lineage partition differs")
if spec["power"]["fibrosis"]["prevalence_is_not_the_null"] is not True:
    raise SystemExit("AUPRC null reference differs")
freeze_tree(root, {
    "artifact_class": "donor_lineage_phenotype_campaign_spec",
    "dataset_id": "gse296875",
    "donors": 39,
    "unit_of_inference": "donor",
    "role": "secondary_development_only",
    "champion_eligible": False,
    "frozen_before_any_fitting": True,
    "outcomes_read": False,
    "model_training_activated": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE296875 phenotype campaign spec: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
