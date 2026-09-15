#!/usr/bin/env bash
#SBATCH --job-name=model-check-010
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
TARGET="${ROOT}/executions/ranking-metric-null-reference-note-v3-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite methods note v3 artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
mkdir -p "${STAGE}/validation"

cat > "${STAGE}/validation/cited_artifacts.sha256" <<'SUMS'
5682a45aef72038e9a8460fbe4b3b43bff14bbb75ca950e3a0292b72a0380b36  executions/ranking-metric-null-reference-note-21107934/ARTIFACTS.json
d4f7e84ab37e5bd9f91e8c67ffe894e6f5184fc16baa2b387b300f96a9ddd328  executions/ranking-metric-null-reference-note-v2-21108960/ARTIFACTS.json
d92205364f18149070b7e636c3cdedc57158feb8513751cfc5c15593bcc0edc7  executions/data-eval-247-21092325/ARTIFACTS.json
f429021ebe8f742314a2073747267047447d07398a3e34eb799b3728c9d829d0  executions/gse274114-activation-21066117/ARTIFACTS.json
92a22bbbf2e8456288c6bc92e2704799da261672d943b033e1de1354183eeebd  executions/model-scoring-078-21092130/ARTIFACTS.json
e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580  executions/gse267145-authoritative-join-21064930/ARTIFACTS.json
9947e411225c3ae3179908fe1cf50bdf9095dd0ac13fda709195d1f85b58b53d  executions/model-check-249-21096291/ARTIFACTS.json
8ae4d2854a9d644141b3d47583feae52465431d1883e4b96b50fe205085e0335  executions/gc-matched-evaluation-windows-21063848/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/cited_artifacts.sha256"

PYTHONPATH=src python3 -m unittest -v \
  tests.unit.test_auprc_reference \
  tests.unit.test_ranking_metric_floor_requirement \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH=src python3 scripts/freeze_ranking_metric_null_reference_note_v3.py \
  --root "${ROOT}" \
  --output "${STAGE}/note" \
  > "${STAGE}/validation/note_stdout.json"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true

sha256sum \
  scripts/freeze_ranking_metric_null_reference_note_v3.py \
  src/masld_bench/evaluators/auprc_reference.py \
  config/campaigns/ranking_metric_floor_requirement_20260825.json \
  config/campaigns/tied_baseline_floor_annotations_20260825.json \
  slurm/freeze_ranking_metric_null_reference_note_v3_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
note = json.loads((root / "note" / "methods_note_v3.json").read_text())
guard = note["computing_a_reference_is_not_recomputing_a_metric"]
if guard["candidate_metrics_recomputed"] != 0:
    raise SystemExit("no candidate metric may be recomputed here")
if guard["candidate_metrics_restated_or_revised"] != 0:
    raise SystemExit("no candidate metric may be restated or revised here")
if note["scope"] != "benchmark_wide_methods_note_not_cohort_specific":
    raise SystemExit("the note must stay phrased for reuse")
if note["supersedes"]["artifacts_sha256"] != (
    "d4f7e84ab37e5bd9f91e8c67ffe894e6f5184fc16baa2b387b300f96a9ddd328"
):
    raise SystemExit("version two binding differs")
if "correction" not in note["supersedes"]:
    raise SystemExit("version three must state what it corrects")
floors = note["computed_floors"]["gse267145"]
if "the_two_tie_mechanisms_are_not_the_same" not in floors:
    raise SystemExit("the corrected tie claim is missing")
if "gain metric" not in note["claim_3_gain_metrics_are_preferable_at_small_n"]["statement"]:
    raise SystemExit("the gain-metric claim is missing")
freeze_tree(root, {
    "artifact_class": "benchmark_methods_note",
    "topic": "ranking_metric_null_references",
    "version": 3,
    "scope": "benchmark_wide",
    "supersedes": "executions/ranking-metric-null-reference-note-21107934",
    "candidate_metrics_recomputed": 0,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Ranking-metric null reference note v2: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
