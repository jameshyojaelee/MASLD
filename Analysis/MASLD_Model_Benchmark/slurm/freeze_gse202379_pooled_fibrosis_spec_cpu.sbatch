#!/usr/bin/env bash
#SBATCH --job-name=model-check-015
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/gse202379-pooled-fibrosis-spec-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"; FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite spec artifact: %s\n' "${TARGET}" >&2; exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
mkdir -p "${STAGE}/validation"
cat > "${STAGE}/validation/bound_fixtures.sha256" <<'SUMS'
6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140  executions/gse296875-phenotype-endpoints-20260825/ARTIFACTS.json
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
3f371ba28847f28c5c17295d50dc6ca65f51c4bde6039ee31e2bcf764603d40f  executions/ranking-metric-null-reference-note-v3-21109019/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

PYTHONPATH=src python3 -m unittest -v \
  tests.unit.test_detectable_effect tests.unit.test_auprc_reference \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH=src python3 scripts/freeze_gse202379_pooled_fibrosis_spec.py \
  --root "${ROOT}" --output "${STAGE}/spec" > "${STAGE}/validation/spec_stdout.json"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum scripts/freeze_gse202379_pooled_fibrosis_spec.py \
          src/masld_bench/evaluators/detectable_effect.py \
          slurm/freeze_gse202379_pooled_fibrosis_spec_cpu.sbatch > "${STAGE}/source.sha256"

PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
spec = json.loads((root / "spec" / "pooled_fibrosis_spec.json").read_text())
g = spec["graded_set"]
if g["graded_donors"] != 40 or g["donors"] != 47 or g["ungraded_donors"] != 7:
    raise SystemExit("graded census differs from the primary source")
if g["distribution"] != {"F0": 3, "F1": 9, "F2": 12, "F3": 12, "F4": 4}:
    raise SystemExit("SAF distribution differs")
if spec["frozen_before_any_gse202379_molecular_data_was_opened"] is not True:
    raise SystemExit("spec must be frozen before molecular data is opened")
if spec["combination_rule"]["naive_pooling_prohibited"] is not True:
    raise SystemExit("naive pooling must be prohibited")
if spec["harmonisation"]["forty_six_is_never_used"] is not True:
    raise SystemExit("the 46-donor table set must be excluded")
if spec["power"]["derived_from_the_frozen_design_not_from_any_result"] is not True:
    raise SystemExit("power must be design-derived")
freeze_tree(root, {
    "artifact_class": "gse202379_pooled_fibrosis_prespecification",
    "datasets": ["gse202379", "gse296875"],
    "unit_of_inference": "donor",
    "graded_donors": 40,
    "frozen_before_molecular_data_opened": True,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE202379 pooled fibrosis spec: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
