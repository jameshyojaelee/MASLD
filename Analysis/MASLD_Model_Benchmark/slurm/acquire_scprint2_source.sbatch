#!/usr/bin/env bash
#SBATCH --job-name=model-source-302
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
commit="7ee2de5aaa3492f326f540f703d44ddeaf1f2fd4"
repository="https://github.com/cantinilab/scPRINT-2.git"
target="${root}/executions/scprint2-source-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite scPRINT-2 source acquisition: %s\n' "${target}" >&2
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

git --git-dir="${stage}/repository.git" archive \
  --format=tar \
  --prefix=source/ \
  --output="${stage}/source.tar" \
  "${commit}"
tar -xf "${stage}/source.tar" -C "${stage}"

cat > "${stage}/expected_source_records.sha256" <<'EOF'
ae49d34bb24814b78a8c3390342986b3df683fe719864f9d6c4569ad27f0ec9a  source/LICENSE
42cf19d5b74199bef4328838910f168f2b00aa0c3cc1fd285d517860a710cfb1  source/README.md
9a24613df23178ff23694fede78172047a116449b54d84103c4baf411633b83a  source/pyproject.toml
62d93b7015491e9534d8768470c46b3b2ae3f9b2ab893079b91dbb84ef546118  source/config/base_v2.yml
fe60cfbf31e7d322859fb01140daecc9ddf23005e201cdbe27cfc55a9be16840  source/data/README.md
b08acbdafa16f39e006ecf5441dd0604a0d8938b7bfb35c1618a7b0974e76e66  source/scprint2/tasks/cell_emb.py
7866483f99ab50886c12a18bfedfff57bfef0993fe80b28d706136005f73103a  source/scprint2/model/model.py
EOF
(
  cd "${stage}"
  sha256sum --check expected_source_records.sha256
  find source -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum \
    > source_tree.sha256
)

source_tree_ledger_sha256="$(sha256sum "${stage}/source_tree.sha256" | cut -d' ' -f1)"
source_tar_sha256="$(sha256sum "${stage}/source.tar" | cut -d' ' -f1)"
license_sha256="$(sha256sum "${stage}/source/LICENSE" | cut -d' ' -f1)"
python3 - "${stage}/ARTIFACTS.json" "${commit}" "${source_tar_sha256}" \
  "${source_tree_ledger_sha256}" "${license_sha256}" <<'PY'
import json
from pathlib import Path
import sys

output, commit, source_tar, source_tree, license_sha = sys.argv[1:]
payload = {
    "checkpoint_deserialized": False,
    "commit": commit,
    "commit_verified": True,
    "license": "GPL-3.0-or-later",
    "license_sha256": license_sha,
    "repository": "cantinilab/scPRINT-2",
    "schema_version": "masld-bench-pinned-source-acquisition-v1",
    "source_release": "1.0.3",
    "source_tar_sha256": source_tar,
    "source_tree_ledger_sha256": source_tree,
    "status": "pass",
}
Path(output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
(
  cd "${stage}"
  sha256sum ARTIFACTS.json expected_source_records.sha256 source.tar source_tree.sha256 \
    > SHA256SUMS
)
chmod -R a-w "${stage}"
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
