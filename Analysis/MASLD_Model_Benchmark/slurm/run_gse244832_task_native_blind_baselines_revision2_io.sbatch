#!/usr/bin/env bash
#SBATCH --job-name=model-train-341
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG=${ROOT}/config/campaigns/gse244832_task_native_blind_baselines_revision2_20260825.json
FIT_PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SOURCE=${ROOT}/executions/model-data-341-${SLURM_JOB_ID}
SOURCE_STAGE=${SOURCE}.staging
RUN=${ROOT}/executions/model-train-341-${SLURM_JOB_ID}
RUN_STAGE=${RUN}.staging

if [[ -e "${SOURCE}" || -e "${SOURCE_STAGE}" || -e "${RUN}" || -e "${RUN_STAGE}" ]]; then
  printf 'Refusing to overwrite a baseline campaign artifact.\n' >&2
  exit 1
fi

cd "${ROOT}"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

sha256sum -c <<'EOF'
ded2530b6ae69eeb7168a347082421c67a3c2a08498498f627d1bc202de5729e  config/campaigns/gse244832_task_native_blind_baselines_revision2_20260825.json
ad68dadfbd84963a6273d3f3312fb2f2ab45e269eacbda9d48e87f2e61371dc9  src/masld_bench/gse244832_task_native_baselines.py
7950821600d6ad3736e458c709044f3cfd391f2f42d15f4681b141179305c5bd  scripts/materialize_gse296875_gse244832_source_window_counts.py
e8b9c0224222898e033640c66354207ecd38a3a40102ab1b4a6f6262a497b1a9  scripts/run_gse244832_task_native_blind_baselines.py
eb536757c0ca9121acf11605082183ddb37060ac9599aa6a02d3e92e8e0a4fcc  scripts/verify_gse244832_task_native_blind_baselines.py
29f913b8736012462dddc92d72304354ee9c4c92b1ec0941a601ef097e09bbee  scripts/commit_gse244832_atac_transport_predictions.py
60dc370af73de685de5dbbf77c7b3749e533030c425cfa35d6c2030e13bdf16b  scripts/py310_compat/tomllib.py
1cbcde8af1a8cd530ff04d93ac55822916cf958d9b3b7f990ca5b7237580c806  tests/unit/test_gse244832_task_native_baselines.py
51fd2323f16b0b951246b04acbd67f51fffceef4fd4fc1ad0bc4d9ecfd525413  tests/unit/test_materialize_gse296875_gse244832_source_window_counts.py
4e31dc1a5730e2b10b4269a2afddb0d827128a83218e5cb69a9878683f74f4c2  tests/unit/test_run_gse244832_task_native_blind_baselines.py
c368d53dd934f9d75050edba7dfe3a35706e64202cdda9a4e868710c0f7538e2  tests/unit/test_verify_gse244832_task_native_blind_baselines.py
194a0a9d96ce036ec1e67409dfd59623c079b706ca3547accf67d9f9278c0a6f  executions/model-check-337-21101061/ARTIFACTS.json
ff1a506d5cca81382911355978097ecd49fcb1b2b48826e8e760a5f2ad1ad35e  executions/gse244832-reference-guarded-atac-exchange-axis-20260825/ARTIFACTS.json
801411d062f9d9f199c0468b6c8c6256263af85ea2805d75e3d166299d190742  executions/gse244832-label-free-atac-query-valid-context-20260825/ARTIFACTS.json
970992d3dc84f2d2c6c2c5357149579fa09b16051ddafb7eed0c1f2949a9ee8a  executions/gse244832-label-free-atac-query-test-context-20260825/ARTIFACTS.json
01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f  executions/gse296875-deduplicated-tn5-bigwigs-21063992/ARTIFACTS.json
ca8f157669d5a6953765525b62ab593db1d4a352e99a73ed1b6db353aab35ef4  executions/gse296875-donor-lineage-fragments-21063889/ARTIFACTS.json
EOF

PYTHONPATH=src:. "${FIT_PYTHON}" -m unittest \
  tests.unit.test_gse244832_task_native_baselines \
  tests.unit.test_materialize_gse296875_gse244832_source_window_counts \
  tests.unit.test_run_gse244832_task_native_blind_baselines \
  tests.unit.test_verify_gse244832_task_native_blind_baselines \
  tests.unit.test_commit_gse244832_atac_transport_predictions

module purge
module load pyBigWig/0.3.18-foss-2022a
PYTHONPATH=scripts/py310_compat:src:. python - <<'PY'
import pyBigWig
import tomllib
from masld_bench.artifacts import freeze_tree

assert pyBigWig.__version__ == "0.3.18"
assert tomllib.loads("value = 1")["value"] == 1
assert callable(freeze_tree)
print("Python 3.10 pyBigWig compatibility bridge passed")
PY
PYTHONPATH=scripts/py310_compat:src:. python scripts/materialize_gse296875_gse244832_source_window_counts.py \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${SOURCE_STAGE}" \
  --workers "${SLURM_CPUS_PER_TASK}"
mv "${SOURCE_STAGE}" "${SOURCE}"
SOURCE_SHA=$(sha256sum "${SOURCE}/ARTIFACTS.json" | cut -d' ' -f1)

module purge
PYTHONPATH=src:. "${FIT_PYTHON}" scripts/run_gse244832_task_native_blind_baselines.py \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --source-input "${SOURCE}" \
  --source-artifacts-sha256 "${SOURCE_SHA}" \
  --output "${RUN_STAGE}"
mv "${RUN_STAGE}" "${RUN}"

sha256sum "${SOURCE}/ARTIFACTS.json" "${RUN}/ARTIFACTS.json"
