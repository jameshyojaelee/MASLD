#!/usr/bin/env bash
#SBATCH --job-name=scprint-preflight
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
runtime="/gpfs/commons/home/jameslee/micromamba/envs/scprint310"
checkpoint="${root}/executions/scprint-acquisition-21015746/v2-medium.ckpt"
target="${root}/executions/scprint-v1-runtime-preflight-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite scPRINT runtime preflight: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export LD_LIBRARY_PATH="${runtime}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

"${runtime}/bin/python" "${root}/scripts/preflight_scprint_v1_runtime.py" \
  --checkpoint "${checkpoint}" \
  --output "${stage}/runtime_preflight.json"
"${runtime}/bin/python" -m pip freeze --all > "${stage}/pip_freeze.txt"
/gpfs/commons/home/jameslee/.local/bin/micromamba list \
  --prefix "${runtime}" --explicit > "${stage}/conda_explicit.txt"

python_sha256="$(sha256sum "${root}/scripts/preflight_scprint_v1_runtime.py" | cut -d' ' -f1)"
runtime_preflight_sha256="$(sha256sum "${stage}/runtime_preflight.json" | cut -d' ' -f1)"
preflight_status="$("${runtime}/bin/python" -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["status"])' \
  "${stage}/runtime_preflight.json")"
cat > "${stage}/ARTIFACTS.json" <<EOF
{
  "outcomes_or_labels_read": false,
  "preflight_source_sha256": "${python_sha256}",
  "runtime_preflight_sha256": "${runtime_preflight_sha256}",
  "schema_version": "masld-bench-scprint-v1-runtime-preflight-artifact-v1",
  "status": "${preflight_status}"
}
EOF
(
  cd "${stage}"
  sha256sum ARTIFACTS.json conda_explicit.txt pip_freeze.txt runtime_preflight.json > SHA256SUMS
)
chmod -R a-w "${stage}"
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
