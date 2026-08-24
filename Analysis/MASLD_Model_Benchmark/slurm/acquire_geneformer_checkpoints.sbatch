#!/usr/bin/env bash
#SBATCH --job-name=masld-geneformer-acquire
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1

module load python/3.11.5-GCCcore-13.2.0

package_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
repository="https://huggingface.co/ctheodoris/Geneformer/resolve/04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5"
target="${package_root}/executions/geneformer-acquisition-${SLURM_JOB_ID}"
stage="${target}.staging"
source_root="${stage}/source"
bundle_root="${stage}/bundles"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite Geneformer acquisition: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${source_root}" "${bundle_root}"

download() {
  local relative="$1"
  local destination="${source_root}/${relative}"
  install -d -m 0750 "$(dirname "${destination}")"
  if [[ -e "${destination}" ]]; then
    printf 'Refusing to overwrite downloaded artifact: %s\n' "${destination}" >&2
    exit 1
  fi
  curl \
    --fail \
    --location \
    --retry 5 \
    --retry-all-errors \
    --show-error \
    --silent \
    --output "${destination}" \
    "${repository}/${relative}?download=true"
}

artifacts=(
  "Geneformer-V1-10M/config.json"
  "Geneformer-V1-10M/model.safetensors"
  "Geneformer-V2-104M/config.json"
  "Geneformer-V2-104M/generation_config.json"
  "Geneformer-V2-104M/model.safetensors"
  "Geneformer-V2-316M/config.json"
  "Geneformer-V2-316M/generation_config.json"
  "Geneformer-V2-316M/model.safetensors"
  "geneformer/ensembl_mapping_dict_gc104M.pkl"
  "geneformer/gene_dictionaries_30m/ensembl_mapping_dict_gc30M.pkl"
  "geneformer/gene_dictionaries_30m/gene_median_dictionary_gc30M.pkl"
  "geneformer/gene_dictionaries_30m/gene_name_id_dict_gc30M.pkl"
  "geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
  "geneformer/gene_median_dictionary_gc104M.pkl"
  "geneformer/gene_name_id_dict_gc104M.pkl"
  "geneformer/token_dictionary_gc104M.pkl"
)
for artifact in "${artifacts[@]}"; do
  download "${artifact}"
done

expected="${stage}/EXPECTED_SHA256SUMS"
cat > "${expected}" <<'EOF'
9cf69ca3bdb0215c4188b54c451b6f02adfe68b8f66011a57d0f32845133fd4b  Geneformer-V1-10M/config.json
a5e33a757431643b3697de7ef6127950cdc49e06e58d4266b3a3ab191b683f14  Geneformer-V1-10M/model.safetensors
467d4492f0dd53b4d60afffe20812db484ca1cf9fdbeb6a6e060e93564f70859  Geneformer-V2-104M/config.json
bf1ba320993168b0d6ee5ff170df2108454d239b28eaac5948d30ee6993cccc1  Geneformer-V2-104M/generation_config.json
fff5cba29ddd8792991fa77b4872246fbe548a178cebda3775cdc72b67780e7f  Geneformer-V2-104M/model.safetensors
2cc4af3442644e84af71814a607b61c958d7b4e5ebde6791ab77b5f534ac6f6e  Geneformer-V2-316M/config.json
bf1ba320993168b0d6ee5ff170df2108454d239b28eaac5948d30ee6993cccc1  Geneformer-V2-316M/generation_config.json
965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3  Geneformer-V2-316M/model.safetensors
0819bcbd869cfa14279449b037eb9ed1d09a91310e77bd1a19d927465030e95c  geneformer/ensembl_mapping_dict_gc104M.pkl
eac0fb0b3007267871b6305ac0003ceba19d4f28d85686cb9067ecf142787869  geneformer/gene_dictionaries_30m/ensembl_mapping_dict_gc30M.pkl
b3b589bb5ec75040d05fc44dd6bf0184cf87f3c362cf158d196a6ed3b7fe5f39  geneformer/gene_dictionaries_30m/gene_median_dictionary_gc30M.pkl
55e67962e79c0039a6c32d43c5c99f38e51964bbcfa32f736150ee1e285c438c  geneformer/gene_dictionaries_30m/gene_name_id_dict_gc30M.pkl
ab9dc40973fa5224d77b793e2fd114cacf3d08423ed9c4c49caf0ba9c7f218f1  geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl
a51c53f6a771d64508dfaf61529df70e394c53bd20856926117ae5d641a24bf5  geneformer/gene_median_dictionary_gc104M.pkl
fabfa0c2f49c598c59ae432a32c3499a5908c033756c663b5e0cddf58deea8e1  geneformer/gene_name_id_dict_gc104M.pkl
67c445f4385127adfc48dcc072320cd65d6822829bf27dd38070e6e787bc597f  geneformer/token_dictionary_gc104M.pkl
EOF

(
  cd "${source_root}"
  sha256sum -c "${expected}"
)

cd "${package_root}"
for model_id in \
  geneformer_v1_10m \
  geneformer_v2_104m \
  geneformer_v2_316m; do
  PYTHONPATH=src python -m masld_bench.cli checkpoint stage-geneformer \
    --model-id "${model_id}" \
    --source-root "${source_root}" \
    --output-root "${bundle_root}"
done

cp "${expected}" "${source_root}/SHA256SUMS"
chmod 0440 "${source_root}/SHA256SUMS"
find "${source_root}" -type f -exec chmod 0440 {} +
find "${bundle_root}" -type f -exec chmod 0440 {} +
mv "${stage}" "${target}"
printf 'Geneformer acquisition candidate: %s\n' "${target}"
