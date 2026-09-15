#!/usr/bin/env bash
#SBATCH --job-name=model-check-1207
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
FIXTURE=${ROOT}/executions/model-data-064-21079902/fixture
ACTIVATION=${ROOT}/executions/model-data-061-21079623/activation
FITVIEWS=${ROOT}/executions/model-training-067-21080776/fit_views
SCORING=${ROOT}/executions/model-scoring-078-21092130

TAG=${TAG:?TAG must be set}
PERMS=${PERMS:-200}
DRAWS=${DRAWS:-50}
BOOTS=${BOOTS:-10000}

TARGET=${ROOT}/executions/gse267145_chromatin_severity_ordering-${TAG}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

"${PYTHON}" scripts/test_gse267145_chromatin_recovers_rna_severity_ordering.py \
  --fixture "${FIXTURE}" \
  --activation "${ACTIVATION}" \
  --fit-views "${FITVIEWS}" \
  --scoring "${SCORING}" \
  --permutations "${PERMS}" \
  --random-draws "${DRAWS}" \
  --bootstrap-replicates "${BOOTS}" \
  --output "${STAGE}/analysis" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/test_gse267145_chromatin_recovers_rna_severity_ordering.py \
  scripts/fit_gse267145_histology_baselines.py \
  slurm/test_gse267145_chromatin_severity_ordering_cpu.sbatch \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
report = json.loads((root / "analysis" / "report.json").read_text())
freeze_tree(root / "analysis", {
    "artifact_class": "gse267145_chromatin_vs_rna_severity_ordering",
    "samples": report["samples"],
    "claim_grade": report["claim_grade"],
    "champion_eligible": False,
    "status": "passed",
})
verify_frozen_tree(root / "analysis")
freeze_tree(root, {
    "artifact_class": "gse267145_chromatin_vs_rna_severity_ordering_campaign",
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'TARGET=%s\n' "${TARGET}"
