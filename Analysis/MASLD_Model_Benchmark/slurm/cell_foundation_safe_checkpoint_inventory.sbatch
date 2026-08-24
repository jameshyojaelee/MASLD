#!/usr/bin/env bash
#SBATCH --job-name=masld-cellfm-inventory
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
model_id="${1:?usage: sbatch cell_foundation_safe_checkpoint_inventory.sbatch MODEL_ID}"

case "${model_id}" in
  regformer)
    checkpoint="${root}/executions/regformer-acquisition-21066043/best_model.pt"
    expected_sha256="8e8e751a0b05caa4e6e4d6b5c0f8ee847ae607757f0431b4bee90cd9ef26a434"
    expected_size="197677201"
    ;;
  scprint2_small_v2)
    checkpoint="${root}/executions/scprint2-acquisition-21066044/18hebyht-final-small.ckpt"
    expected_sha256="2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822"
    expected_size="889706878"
    ;;
  scgpt_whole_human)
    checkpoint="${root}/executions/scgpt-acquisition-21015743/scgpt_whole_human/best_model.pt"
    expected_sha256="6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922"
    expected_size="205385258"
    ;;
  scgpt_continual)
    checkpoint="${root}/executions/scgpt-acquisition-21015743/scgpt_continual/best_model.pt"
    expected_sha256="ad0252a1971e0cd619b7116dbab3177432236c4537225d54280a2aa7e5fe402a"
    expected_size="207861754"
    ;;
  uce_4l)
    checkpoint="${root}/executions/uce-acquisition-21015745/4layer_model.torch"
    expected_sha256="acb28f3f0a1d803e4a4ffe891b9bab38bf93c84762dc06b2452f0d515da91560"
    expected_size="3403514339"
    ;;
  uce_33l)
    checkpoint="${root}/executions/uce-acquisition-21015745/33l_8ep_1024t_1280.torch"
    expected_sha256="3f458726196308e171611ed28394b55865708749f82af68cfd7771d5dfee661e"
    expected_size="5686283829"
    ;;
  uce_all_tokens)
    checkpoint="${root}/executions/uce-acquisition-21015745/all_tokens.torch"
    expected_sha256="07397ab3828502fb7d0bab658125c47f145f989f384db8f55fe106f6826b2a54"
    expected_size="2979205876"
    ;;
  scimilarity_v1_1_encoder)
    checkpoint="${root}/executions/scimilarity-model-20260823/model_v1.1/encoder.ckpt"
    expected_sha256="a08f023788bdbded191e2d9344e5d4240ffc7d7505eca4e0b7a37526b858d19f"
    expected_size="124616535"
    ;;
  scprint_v1_5_medium)
    checkpoint="${root}/executions/scprint-acquisition-21015746/v2-medium.ckpt"
    expected_sha256="a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5"
    expected_size="221775592"
    ;;
  *)
    printf 'Unsupported model identity: %s\n' "${model_id}" >&2
    exit 2
    ;;
esac

target="${root}/executions/${model_id}-safe-inventory-${SLURM_JOB_ID}"
stage="${target}.staging"
if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite inventory candidate: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

python3 "${root}/scripts/torch_checkpoint_safe_inventory.py" \
  --input "${checkpoint}" \
  --output "${stage}/checkpoint_inventory.json" \
  --expected-sha256 "${expected_sha256}" \
  --expected-size "${expected_size}"

python3 --version > "${stage}/python_version.txt" 2>&1
python3 -m pip freeze > "${stage}/pip_freeze.txt"
(
  cd "${stage}"
  sha256sum checkpoint_inventory.json python_version.txt pip_freeze.txt > SHA256SUMS
)
sha256sum \
  "${root}/scripts/torch_checkpoint_safe_inventory.py" \
  "${root}/slurm/cell_foundation_safe_checkpoint_inventory.sbatch" \
  > "${stage}/source_manifest.sha256"
chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'Safe checkpoint inventory candidate: %s\n' "${target}"
