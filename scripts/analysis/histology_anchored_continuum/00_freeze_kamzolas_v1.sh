#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${HAC_OUT_ROOT:-}" ]]; then
  echo "HAC_OUT_ROOT is unset" >&2
  exit 1
fi

repo_url="https://github.com/kamzolas/MASLD---Continuous-trajectory-approach.git"
expected_commit="cbcf7764b625d851ca4e2e58b400b921cfb8f39a"
external_root="${HAC_OUT_ROOT}/external"
checkout="${external_root}/kamzolas_v1"

if [[ -e "${checkout}" ]]; then
  echo "Refusing to overwrite ${checkout}" >&2
  exit 1
fi

mkdir -p "${external_root}"
git clone --branch v1.0 --depth 1 "${repo_url}" "${checkout}"
observed_commit="$(git -C "${checkout}" rev-parse HEAD)"
if [[ "${observed_commit}" != "${expected_commit}" ]]; then
  echo "Kamzolas v1 commit mismatch: expected ${expected_commit}; observed ${observed_commit}" >&2
  exit 1
fi
if [[ -n "$(git -C "${checkout}" status --porcelain)" ]]; then
  echo "Frozen checkout is dirty" >&2
  exit 1
fi

printf '%s\t%s\t%s\n' "tag" "commit" "repository" > "${external_root}/kamzolas_release.tsv"
printf '%s\t%s\t%s\n' "v1.0" "${observed_commit}" "${repo_url}" >> "${external_root}/kamzolas_release.tsv"
printf '%s\n' "${checkout}"
