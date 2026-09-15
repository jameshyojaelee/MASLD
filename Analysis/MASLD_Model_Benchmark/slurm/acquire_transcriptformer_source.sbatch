#!/usr/bin/env bash
#SBATCH --job-name=tf-source
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
commit="c943a89a4de9511a8ab1010715bf1cfa8827cc99"
repository="https://github.com/czi-ai/transcriptformer.git"
target="${root}/executions/transcriptformer-source-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite source acquisition: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

git init --bare "${stage}/repository.git"
git --git-dir="${stage}/repository.git" remote add origin "${repository}"
git --git-dir="${stage}/repository.git" fetch --depth=1 origin "${commit}"
observed_commit="$(git --git-dir="${stage}/repository.git" rev-parse FETCH_HEAD^{commit})"
if [[ "${observed_commit}" != "${commit}" ]]; then
  printf 'Commit mismatch: expected %s, observed %s\n' \
    "${commit}" "${observed_commit}" >&2
  exit 1
fi
git ls-remote --tags "${repository}" \
  refs/tags/v0.6.0 'refs/tags/v0.6.0^{}' > "${stage}/release_tag_refs.tsv"

git --git-dir="${stage}/repository.git" archive \
  --format=tar \
  --prefix=source/ \
  --output="${stage}/source.tar" \
  "${commit}"
tar -xf "${stage}/source.tar" -C "${stage}"

cat > "${stage}/expected_source_records.sha256" <<'EOF'
aaf473b176a223eeccd843bad7d5453501d826f6f939d44f8957b8f6f2e5cc81  source/LICENSE.md
15d1a28e0af8565e2337ef8fda144444b9c0518f1fb4f8bb66e3da8aac7431f0  source/README.md
20c8bc58a9069b2835b08a5c35ce906740aad923519aeb9e1e97f75689a96583  source/pyproject.toml
45f51bbcb76fcc6a1454fd8f249c814070b05020bee16ab9dc6f823e24ed4c29  source/src/transcriptformer/cli/download_artifacts.py
4e8b3b04e3ede56c548e3eacfcb5ece22232835606424fa3efc21f78bf7eccb7  source/src/transcriptformer/data/dataloader.py
45ea496cfb7203eb6f7f4f262f9502e5b2986b7349b7261d4bef67ca870dda0a  source/conf/inference_config.yaml
EOF
(
  cd "${stage}"
  sha256sum --check expected_source_records.sha256
  find source -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum \
    > source_tree.sha256
)

source_tree_ledger_sha256="$(sha256sum "${stage}/source_tree.sha256" | cut -d' ' -f1)"
source_tar_sha256="$(sha256sum "${stage}/source.tar" | cut -d' ' -f1)"
license_sha256="$(sha256sum "${stage}/source/LICENSE.md" | cut -d' ' -f1)"
cat > "${stage}/ARTIFACTS.json" <<EOF
{
  "commit": "${commit}",
  "commit_verified": true,
  "license": "MIT",
  "license_sha256": "${license_sha256}",
  "source_version": "0.6.1",
  "v0_6_0_tag_is_pinned_commit": false,
  "repository": "czi-ai/transcriptformer",
  "schema_version": "masld-bench-pinned-source-acquisition-v1",
  "source_tar_sha256": "${source_tar_sha256}",
  "source_tree_ledger_sha256": "${source_tree_ledger_sha256}",
  "status": "pass"
}
EOF
(
  cd "${stage}"
  sha256sum \
    ARTIFACTS.json \
    expected_source_records.sha256 \
    release_tag_refs.tsv \
    source.tar \
    source_tree.sha256 \
    > SHA256SUMS
)
chmod -R a-w "${stage}"
mv "${stage}" "${target}"
printf 'Pinned TranscriptFormer source candidate: %s\n' "${target}"
