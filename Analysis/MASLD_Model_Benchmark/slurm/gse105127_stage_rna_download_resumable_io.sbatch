#!/usr/bin/env bash
#SBATCH --job-name=data-stage-220
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-7
#SBATCH --output=executions/logs/data-stage-220-%A_%a.out
#SBATCH --error=executions/logs/data-stage-220-%A_%a.err

source slurm/gse105127_production_common.sh
RESUMABLE_SCRIPT=scripts/download_gse105127_rna_bundle_resumable.py
RESUMABLE_SCRIPT_SHA256=2b93452ecbefbb38a64721b0a736b90fd1f126e5bc5d04b6c977fd563c5416fe
VALIDATION=${ROOT}/executions/data-check-221-21087852
VALIDATION_SHA256=e1d2ebb5096c8ac3ba2ba7b4176edae77325a04c0ba949e9664e423013170617
printf '%s  %s\n' "${RESUMABLE_SCRIPT_SHA256}" "${RESUMABLE_SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' "${VALIDATION_SHA256}" "${VALIDATION}/ARTIFACTS.json" | sha256sum --check --strict
"${PYTHON}" - "${VALIDATION}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
"${PYTHON}" -m unittest -v tests.unit.test_download_gse105127_rna_bundle_resumable
begin_attempt rna-download-resumable
"${PYTHON}" "${RESUMABLE_SCRIPT}" \
  --plan-root "${PLAN}" \
  --remote-identity-root "${ROOT}/executions/data-check-218-21087748" \
  --bundle-id "${SLURM_ARRAY_TASK_ID}" \
  --output "${RNA_RAW}" \
  --attempt "${ATTEMPT}/members"
finish_attempt
