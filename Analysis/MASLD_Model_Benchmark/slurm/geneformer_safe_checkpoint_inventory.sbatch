#!/usr/bin/env bash
#SBATCH --job-name=masld-gf-inventory
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
source_root="${root}/executions/geneformer-acquisition-21014827/source"
model_id="${1:?usage: sbatch geneformer_safe_checkpoint_inventory.sbatch MODEL_ID}"

case "${model_id}" in
  geneformer_v1_10m)
    checkpoint="${source_root}/Geneformer-V1-10M/model.safetensors"
    expected_sha256="a5e33a757431643b3697de7ef6127950cdc49e06e58d4266b3a3ab191b683f14"
    expected_size="41183536"
    ;;
  geneformer_v2_104m)
    checkpoint="${source_root}/Geneformer-V2-104M/model.safetensors"
    expected_sha256="fff5cba29ddd8792991fa77b4872246fbe548a178cebda3775cdc72b67780e7f"
    expected_size="417571156"
    ;;
  geneformer_v2_316m)
    checkpoint="${source_root}/Geneformer-V2-316M/model.safetensors"
    expected_sha256="965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3"
    expected_size="1265455076"
    ;;
  *)
    printf 'Unsupported Geneformer model identity: %s\n' "${model_id}" >&2
    exit 2
    ;;
esac

target="${root}/executions/${model_id}-safe-inventory-${SLURM_JOB_ID}"
stage="${target}.staging"
if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite inventory candidate: %s\n' "${target}" >&2
  exit 1
fi
(
  cd "${source_root}"
  sha256sum -c SHA256SUMS
)
install -d -m 0750 "${stage}"

python3 "${root}/scripts/safetensors_safe_inventory.py" \
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
  "${root}/scripts/safetensors_safe_inventory.py" \
  "${root}/slurm/geneformer_safe_checkpoint_inventory.sbatch" \
  > "${stage}/source_manifest.sha256"
chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'Geneformer safe checkpoint inventory: %s\n' "${target}"
