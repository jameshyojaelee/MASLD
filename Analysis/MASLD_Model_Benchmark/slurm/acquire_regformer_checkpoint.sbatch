#!/usr/bin/env bash
#SBATCH --job-name=masld-regformer-acquire
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

package_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
target="${package_root}/executions/regformer-acquisition-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite RegFormer acquisition: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

fetch() {
  local file_id="$1"
  local output="$2"
  curl \
    --fail \
    --location \
    --retry 5 \
    --retry-all-errors \
    --show-error \
    --silent \
    --output "${stage}/${output}.part" \
    "https://ndownloader.figshare.com/files/${file_id}"
  mv "${stage}/${output}.part" "${stage}/${output}"
}

fetch 57423100 best_model.pt
fetch 57423097 vocab.json
fetch 62134482 output_graph_no_cycles.dgl
fetch 62134521 pretrain.toml

observed_bytes="$(stat --format='%s' "${stage}/best_model.pt")"
if [[ "${observed_bytes}" != "197677201" ]]; then
  printf 'Checkpoint byte-size mismatch: expected 197677201, observed %s\n' \
    "${observed_bytes}" >&2
  exit 1
fi
printf '%s  %s\n' \
  "6f0dea8ae1508077d1851cdcda6fce42" "${stage}/best_model.pt" | md5sum -c -
printf '%s  %s\n' \
  "acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f" \
  "${stage}/vocab.json" | sha256sum -c -
printf '%s  %s\n' \
  "f0a530d760c82ef9e264eec3f3f5a230f295fd3c3e3a6b347e2f84a784807397" \
  "${stage}/output_graph_no_cycles.dgl" | sha256sum -c -
printf '%s  %s\n' \
  "f599274fed2ff2d80a96d2e318e40941ad4c67adb08946b72a990e37b25c57ea" \
  "${stage}/pretrain.toml" | sha256sum -c -

(
  cd "${stage}"
  sha256sum best_model.pt vocab.json output_graph_no_cycles.dgl pretrain.toml \
    > SHA256SUMS
)
chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'RegFormer acquisition candidate: %s\n' "${target}"
printf 'The checkpoint was not deserialized or executed. Restricted inspection remains mandatory.\n'
