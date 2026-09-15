#!/usr/bin/env bash
#SBATCH --job-name=model-check-006
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/ranking-metric-null-reference-note-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite methods note artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
mkdir -p "${STAGE}/validation"

cat > "${STAGE}/validation/bound_fixtures.sha256" <<'SUMS'
6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140  executions/gse296875-phenotype-endpoints-20260825/ARTIFACTS.json
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
0fcf833422f0316644885560b58e0a540b8e98071d7c4701c6801706230c1afd  executions/gse296875-phenotype-scores-21107640/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

PYTHONPATH=src python3 -m unittest -v \
  tests.unit.test_auprc_reference \
  tests.unit.test_gse296875_phenotype_scoring \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH=src python3 scripts/freeze_ranking_metric_null_reference_note.py \
  --root "${ROOT}" \
  --output "${STAGE}/note" \
  > "${STAGE}/validation/note_stdout.json"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true

sha256sum \
  scripts/freeze_ranking_metric_null_reference_note.py \
  src/masld_bench/evaluators/auprc_reference.py \
  tests/unit/test_auprc_reference.py \
  slurm/freeze_ranking_metric_null_reference_note_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
note = json.loads((root / "note" / "methods_note.json").read_text())
if note["scope"] != "benchmark_wide_methods_note_not_cohort_specific":
    raise SystemExit("the note must be phrased for reuse, not for one cohort")
if note["module_is_lane_agnostic"] is not True:
    raise SystemExit("the bound module must be lane agnostic")
example = note["worked_example"]["fibrosis"]
if example["observed_auprc_all_below_the_random_score_mean"] is not True:
    raise SystemExit("worked example no longer demonstrates the mechanism")
if example["prevalence"] >= example["continuous_random_score_reference"]["null_mean"]:
    raise SystemExit("prevalence must sit below the random-score mean")
freeze_tree(root, {
    "artifact_class": "benchmark_methods_note",
    "topic": "ranking_metric_null_references",
    "scope": "benchmark_wide",
    "applies_beyond_gse296875": True,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Ranking-metric null reference note: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
