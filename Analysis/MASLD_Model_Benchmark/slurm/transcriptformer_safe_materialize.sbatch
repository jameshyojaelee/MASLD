#!/usr/bin/env bash
#SBATCH --job-name=masld-tf-materialize
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
archive="${root}/executions/transcriptformer-acquisition-21066078/tf_sapiens.tar.gz"
inventory_root="${root}/executions/transcriptformer-safe-inventory-21066109"
inventory="${inventory_root}/safe_archive_inventory.json"
archive_sha256="d108a7a9ee7fdf38e0dcd4e58593e94304cba927df6ec63e82ea74e2b177e529"
checkpoint_sha256="eff027e7393f92a32cadcbb2e13a917e24b62b71341571bf40f38ff6ec6d3021"
checkpoint_size="1719542526"
target="${root}/executions/transcriptformer-materialized-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer materialization: %s\n' "${target}" >&2
  exit 1
fi
(
  cd "${inventory_root}"
  sha256sum -c SHA256SUMS
)
install -d -m 0750 "${stage}"

python3 "${root}/scripts/safe_tar_materialize.py" \
  --archive "${archive}" \
  --inventory "${inventory}" \
  --output "${stage}/bundle" \
  --receipt "${stage}/materialization_receipt.json" \
  --expected-archive-sha256 "${archive_sha256}"

python3 "${root}/scripts/torch_checkpoint_safe_inventory.py" \
  --input "${stage}/bundle/tf_sapiens/model_weights.pt" \
  --output "${stage}/checkpoint_inventory.json" \
  --expected-sha256 "${checkpoint_sha256}" \
  --expected-size "${checkpoint_size}"

python3 --version > "${stage}/python_version.txt" 2>&1
python3 -m pip freeze > "${stage}/pip_freeze.txt"
(
  cd "${stage}"
  sha256sum materialization_receipt.json checkpoint_inventory.json \
    python_version.txt pip_freeze.txt > SHA256SUMS
  find bundle -type f -print0 | sort -z | xargs -0 sha256sum > bundle.sha256
)
sha256sum \
  "${root}/scripts/safe_tar_materialize.py" \
  "${root}/scripts/torch_checkpoint_safe_inventory.py" \
  "${root}/slurm/transcriptformer_safe_materialize.sbatch" \
  > "${stage}/source_manifest.sha256"
chmod 0440 "${stage}"/*.json "${stage}"/*.txt "${stage}"/*.sha256 "${stage}/SHA256SUMS"
mv "${stage}" "${target}"
printf 'TranscriptFormer safe materialization: %s\n' "${target}"
