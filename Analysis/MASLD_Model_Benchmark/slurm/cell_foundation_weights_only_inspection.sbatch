#!/usr/bin/env bash
#SBATCH --job-name=masld-cellfm-inspect
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONNOUSERSITE=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
python_env="/gpfs/commons/home/jameslee/micromamba/envs/decima"
micromamba_bin="/gpfs/commons/home/jameslee/.local/bin/micromamba"
model_id="${1:?usage: sbatch cell_foundation_weights_only_inspection.sbatch MODEL_ID}"

case "${model_id}" in
  regformer)
    checkpoint="${root}/executions/regformer-acquisition-21066043/best_model.pt"
    inventory="${root}/executions/regformer-safe-inventory-21066069"
    expected_sha256="8e8e751a0b05caa4e6e4d6b5c0f8ee847ae607757f0431b4bee90cd9ef26a434"
    expected_size="197677201"
    ;;
  scprint2_small_v2)
    checkpoint="${root}/executions/scprint2-acquisition-21066044/18hebyht-final-small.ckpt"
    inventory="${root}/executions/scprint2_small_v2-safe-inventory-21066070"
    expected_sha256="2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822"
    expected_size="889706878"
    ;;
  transcriptformer_sapiens)
    checkpoint="${root}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens/model_weights.pt"
    inventory="${root}/executions/transcriptformer-materialized-21066178"
    expected_sha256="eff027e7393f92a32cadcbb2e13a917e24b62b71341571bf40f38ff6ec6d3021"
    expected_size="1719542526"
    ;;
  scgpt_continual)
    checkpoint="${root}/executions/scgpt-acquisition-21015743/scgpt_continual/best_model.pt"
    inventory="${root}/executions/scgpt_continual-safe-inventory-21066202"
    expected_sha256="ad0252a1971e0cd619b7116dbab3177432236c4537225d54280a2aa7e5fe402a"
    expected_size="207861754"
    ;;
  scgpt_whole_human)
    checkpoint="${root}/executions/scgpt-acquisition-21015743/scgpt_whole_human/best_model.pt"
    inventory="${root}/executions/scgpt_whole_human-safe-inventory-21066205"
    expected_sha256="6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922"
    expected_size="205385258"
    ;;
  uce_4l)
    checkpoint="${root}/executions/uce-acquisition-21015745/4layer_model.torch"
    inventory="${root}/executions/uce_4l-safe-inventory-21066201"
    expected_sha256="acb28f3f0a1d803e4a4ffe891b9bab38bf93c84762dc06b2452f0d515da91560"
    expected_size="3403514339"
    ;;
  uce_33l)
    checkpoint="${root}/executions/uce-acquisition-21015745/33l_8ep_1024t_1280.torch"
    inventory="${root}/executions/uce_33l-safe-inventory-21066199"
    expected_sha256="3f458726196308e171611ed28394b55865708749f82af68cfd7771d5dfee661e"
    expected_size="5686283829"
    ;;
  uce_all_tokens)
    checkpoint="${root}/executions/uce-acquisition-21015745/all_tokens.torch"
    inventory="${root}/executions/uce_all_tokens-safe-inventory-21066204"
    expected_sha256="07397ab3828502fb7d0bab658125c47f145f989f384db8f55fe106f6826b2a54"
    expected_size="2979205876"
    ;;
  scimilarity_v1_1_encoder)
    checkpoint="${root}/executions/scimilarity-model-20260823/model_v1.1/encoder.ckpt"
    inventory="${root}/executions/scimilarity_v1_1_encoder-safe-inventory-21066203"
    expected_sha256="a08f023788bdbded191e2d9344e5d4240ffc7d7505eca4e0b7a37526b858d19f"
    expected_size="124616535"
    ;;
  scprint_v1_5_medium)
    checkpoint="${root}/executions/scprint-acquisition-21015746/v2-medium.ckpt"
    inventory="${root}/executions/scprint_v1_5_medium-safe-inventory-21066200"
    expected_sha256="a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5"
    expected_size="221775592"
    ;;
  *)
    printf 'Unsupported model identity: %s\n' "${model_id}" >&2
    exit 2
    ;;
esac

target="${root}/executions/${model_id}-weights-only-inspection-${SLURM_JOB_ID}"
stage="${target}.staging"
if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite inspection candidate: %s\n' "${target}" >&2
  exit 1
fi
(
  cd "${inventory}"
  sha256sum -c SHA256SUMS
)
install -d -m 0750 "${stage}"

"${micromamba_bin}" run --prefix "${python_env}" \
  python "${root}/scripts/torch_checkpoint_weights_only_inspect.py" \
  --checkpoint "${checkpoint}" \
  --output "${stage}/weights_only_inspection.json" \
  --expected-sha256 "${expected_sha256}" \
  --expected-size "${expected_size}"

"${micromamba_bin}" run --prefix "${python_env}" python --version \
  > "${stage}/python_version.txt" 2>&1
"${micromamba_bin}" run --prefix "${python_env}" python -m pip freeze \
  > "${stage}/pip_freeze.txt"
(
  cd "${stage}"
  sha256sum weights_only_inspection.json python_version.txt pip_freeze.txt > SHA256SUMS
)
sha256sum \
  "${root}/scripts/torch_checkpoint_weights_only_inspect.py" \
  "${root}/slurm/cell_foundation_weights_only_inspection.sbatch" \
  > "${stage}/source_manifest.sha256"
chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'Restricted weights-only inspection candidate: %s\n' "${target}"
