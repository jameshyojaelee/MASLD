#!/usr/bin/env bash
#SBATCH --job-name=masld-fragment-membership
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ACQUISITION="${ROOT}/executions/gse296875-fragments-direct-21063797"
MAPPING=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/external-gse296875-metadata-v45-2026-08-17T230000Z/mapping_metadata.tsv
LABELS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/external-gse296875-author-labels-v46-2026-08-17T234000Z/author_labels.tsv
TARGET="${ROOT}/executions/gse296875-fragment-membership-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite fragment membership: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' 038bb09e1056d86bb1119137450fdef856f9129cb08cb8b2fdb27735d869a303 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 5d3d946607f55b140ffa4a44369cf35f6e9e91f0f17276afc0e3e9f008554781 \
  "${MAPPING}" | sha256sum --check --strict
printf '%s  %s\n' 45c077e9beb8607402bd2ed7bce0586596430d401f290ca86978f91e230ea67d \
  "${LABELS}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python scripts/build_gse296875_fragment_membership.py \
  --acquisition "${ACQUISITION}" \
  --acquisition-artifacts-sha256 038bb09e1056d86bb1119137450fdef856f9129cb08cb8b2fdb27735d869a303 \
  --mapping "${MAPPING}" \
  --labels "${LABELS}" \
  --source-path scripts/build_gse296875_fragment_membership.py \
  --source-path slurm/build_gse296875_fragment_membership.sbatch \
  --output "${STAGE}"
mv "${STAGE}" "${TARGET}"
printf 'Fragment membership: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
