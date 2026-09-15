#!/usr/bin/env bash
#SBATCH --job-name=model-check-008
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
TARGET="${ROOT}/executions/ranking-metric-floor-requirement-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite floor-requirement artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
mkdir -p "${STAGE}/validation" "${STAGE}/declarations"

# Every output file the declarations cite must still verify before they are frozen.
cat > "${STAGE}/validation/cited_artifacts.sha256" <<'SUMS'
9947e411225c3ae3179908fe1cf50bdf9095dd0ac13fda709195d1f85b58b53d  executions/model-check-249-21096291/ARTIFACTS.json
8ae4d2854a9d644141b3d47583feae52465431d1883e4b96b50fe205085e0335  executions/gc-matched-evaluation-windows-21063848/ARTIFACTS.json
d92205364f18149070b7e636c3cdedc57158feb8513751cfc5c15593bcc0edc7  executions/data-eval-247-21092325/ARTIFACTS.json
92a22bbbf2e8456288c6bc92e2704799da261672d943b033e1de1354183eeebd  executions/model-scoring-078-21092130/ARTIFACTS.json
29925cd23db16b603186868f230c3cf9aff6d4e0aeedebe5163eab217ad77c83  executions/model-eval-084-21096365/ARTIFACTS.json
f429021ebe8f742314a2073747267047447d07398a3e34eb799b3728c9d829d0  executions/gse274114-activation-21066117/ARTIFACTS.json
5682a45aef72038e9a8460fbe4b3b43bff14bbb75ca950e3a0292b72a0380b36  executions/ranking-metric-null-reference-note-21107934/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/cited_artifacts.sha256"

PYTHONPATH=src python3 -m unittest -v tests.unit.test_ranking_metric_floor_requirement \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

cp config/campaigns/ranking_metric_floor_requirement_20260825.json "${STAGE}/declarations/"
cp config/campaigns/tied_baseline_floor_annotations_20260825.json "${STAGE}/declarations/"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true

sha256sum \
  config/campaigns/ranking_metric_floor_requirement_20260825.json \
  config/campaigns/tied_baseline_floor_annotations_20260825.json \
  tests/unit/test_ranking_metric_floor_requirement.py \
  slurm/freeze_ranking_metric_floor_requirement_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
requirement = json.loads(
    (root / "declarations" / "ranking_metric_floor_requirement_20260825.json").read_text()
)
annotations = json.loads(
    (root / "declarations" / "tied_baseline_floor_annotations_20260825.json").read_text()
)
if requirement["applies_prospectively_only"] is not True:
    raise SystemExit("the floor requirement must be prospective only")
if requirement["no_frozen_metric_is_altered"] is not True:
    raise SystemExit("the floor requirement must alter no frozen metric")
if annotations["record_only"] is not True or annotations["consumed_by_runtime"] is not False:
    raise SystemExit("the baseline annotations must be record-only")
if annotations["reported_metrics_unchanged"] is not True:
    raise SystemExit("the baseline annotations must change no reported metric")
for entry in annotations["annotated_baselines"]:
    if entry["is_a_valid_chance_level_floor"] is not False:
        raise SystemExit("a constant baseline may not be recorded as a chance level")
freeze_tree(root, {
    "artifact_class": "ranking_metric_floor_requirement",
    "scope": "benchmark_wide",
    "applies_prospectively_only": True,
    "record_only_annotations_included": True,
    "candidate_metrics_recomputed": 0,
    "frozen_metrics_altered": 0,
    "outcomes_read": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Ranking-metric floor requirement: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
