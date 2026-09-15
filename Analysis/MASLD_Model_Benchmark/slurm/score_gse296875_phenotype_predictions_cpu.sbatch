#!/usr/bin/env bash
#SBATCH --job-name=model-check-004
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

: "${PREDICTIONS:?PREDICTIONS path is required}"
: "${PREDICTIONS_SHA256:?PREDICTIONS_SHA256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/gse296875-phenotype-scores-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite score artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
mkdir -p "${STAGE}/validation"

cat > "${STAGE}/validation/bound_fixtures.sha256" <<SUMS
6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140  executions/gse296875-phenotype-endpoints-20260825/ARTIFACTS.json
e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416  executions/gse296875-donor-folds-20260825/ARTIFACTS.json
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
${PREDICTIONS_SHA256}  ${PREDICTIONS}/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

PYTHONPATH=src python3 -m unittest -v \
  tests.unit.test_gse296875_phenotype_scoring \
  tests.unit.test_gse296875_phenotype_campaign_spec \
  tests.contract.test_gse296875_phenotype \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH=src python3 scripts/score_gse296875_phenotype_predictions.py \
  --root "${ROOT}" \
  --predictions "${PREDICTIONS}" \
  --predictions-sha256 "${PREDICTIONS_SHA256}" \
  --output "${STAGE}/scores" \
  > "${STAGE}/validation/score_stdout.json"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true

sha256sum \
  scripts/score_gse296875_phenotype_predictions.py \
  src/masld_bench/gse296875_phenotype_scoring.py \
  tests/unit/test_gse296875_phenotype_scoring.py \
  slurm/score_gse296875_phenotype_predictions_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "scores" / "audit.json").read_text())
if audit["confirmatory_family_size"] != 12:
    raise SystemExit("confirmatory family size differs from the frozen spec")
if audit["pooled_out_of_fold_metric_computed_once_per_endpoint"] is not True:
    raise SystemExit("pooled metric assertion missing")
if audit["champion_eligible"] is not False:
    raise SystemExit("this cohort cannot confirm a champion")
freeze_tree(root, {
    "artifact_class": "donor_lineage_phenotype_pooled_scores",
    "dataset_id": "gse296875",
    "donors": 39,
    "unit_of_inference": "donor",
    "campaign_spec_sha256": "474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58",
    "champion_eligible": False,
    "role": "secondary_development_only",
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE296875 phenotype scores: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
