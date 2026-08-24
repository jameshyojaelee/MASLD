#!/usr/bin/env bash
#SBATCH --job-name=masld-fragment-split
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ACQUISITION="${ROOT}/executions/gse296875-fragments-direct-21063797"
MEMBERSHIP="${ROOT}/executions/gse296875-fragment-membership-21063829"
TARGET="${ROOT}/executions/gse296875-donor-lineage-fragments-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite donor-lineage fragments: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' 038bb09e1056d86bb1119137450fdef856f9129cb08cb8b2fdb27735d869a303 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64 \
  "${MEMBERSHIP}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python scripts/split_gse296875_donor_lineage_fragments.py \
  --acquisition "${ACQUISITION}" \
  --acquisition-artifacts-sha256 038bb09e1056d86bb1119137450fdef856f9129cb08cb8b2fdb27735d869a303 \
  --membership "${MEMBERSHIP}" \
  --membership-artifacts-sha256 582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64 \
  --split-workers 8 \
  --merge-workers 16 \
  --source-path scripts/split_gse296875_donor_lineage_fragments.py \
  --source-path slurm/split_gse296875_donor_lineage_fragments.sbatch \
  --output "${STAGE}"
mv "${STAGE}" "${TARGET}"
printf 'Donor-lineage fragments: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
