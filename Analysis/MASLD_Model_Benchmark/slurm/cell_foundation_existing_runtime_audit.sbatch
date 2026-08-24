#!/usr/bin/env bash
#SBATCH --job-name=masld-cellfm-rtaudit
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
micromamba_bin="/gpfs/commons/home/jameslee/.local/bin/micromamba"
lane="${1:?usage: sbatch cell_foundation_existing_runtime_audit.sbatch LANE}"

case "${lane}" in
  decima|scimilarity|scprint310|scprint2x)
    python_env="/gpfs/commons/home/jameslee/micromamba/envs/${lane}"
    ;;
  *)
    printf 'Unsupported runtime lane: %s\n' "${lane}" >&2
    exit 2
    ;;
esac
if [[ ! -x "${python_env}/bin/python" ]]; then
  printf 'Existing runtime is absent: %s\n' "${python_env}" >&2
  exit 1
fi

target="${root}/executions/${lane}-existing-cellfm-runtime-audit-${SLURM_JOB_ID}"
stage="${target}.staging"
if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite runtime audit: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

"${micromamba_bin}" list --prefix "${python_env}" --explicit > "${stage}/conda_explicit.txt"
"${micromamba_bin}" run --prefix "${python_env}" python -m pip freeze > "${stage}/pip_freeze.txt"
"${micromamba_bin}" run --prefix "${python_env}" \
  python "${root}/scripts/audit_existing_cell_foundation_runtime.py" \
  --lane "${lane}" --output "${stage}/runtime_audit.json"
nvidia-smi --query-gpu=name,uuid,driver_version,memory.total --format=csv,noheader \
  > "${stage}/nvidia_smi.csv"

sha256sum \
  "${root}/scripts/audit_existing_cell_foundation_runtime.py" \
  "${root}/slurm/cell_foundation_existing_runtime_audit.sbatch" \
  > "${stage}/source_manifest.sha256"
(
  cd "${stage}"
  sha256sum conda_explicit.txt pip_freeze.txt runtime_audit.json nvidia_smi.csv \
    > SHA256SUMS
)
chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'Existing runtime audit: %s\n' "${target}"
