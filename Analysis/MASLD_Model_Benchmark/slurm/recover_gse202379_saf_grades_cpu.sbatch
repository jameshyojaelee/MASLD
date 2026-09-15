#!/usr/bin/env bash
#SBATCH --job-name=SAFrecover
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/logs/saf_recover_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/logs/saf_recover_%j.err

# Additive recovery of GSE202379 SAF (steatosis/activity/fibrosis) donor grades.
# Wall time is 1h rather than the 90h project default: this is a single-pass
# parse of one 200KB gzipped GEO series matrix and completes in under a second.

set -eo pipefail
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
export R_LIBS_USER=""
export PYTHONNOUSERSITE=1
export PYTHONWARNINGS=ignore
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

BENCH="${MASLD_PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark"
EXEC="${BENCH}/executions/gse202379-saf-recovery-${SLURM_JOB_ID}"
LOCK="${EXEC}/saf_lock"
VAL="${EXEC}/validation"
mkdir -p "${LOCK}" "${VAL}"

cd "${BENCH}"
echo "=== SAF recovery | job ${SLURM_JOB_ID} | $(hostname) | $(date) ==="

# 1. recover
python scripts/recover_gse202379_saf_grades.py --outdir "${LOCK}" \
  > "${VAL}/recovery_stdout.json" 2> "${VAL}/recovery_stderr.txt"

# 2. environment
pip freeze > "${VAL}/pip_freeze.txt"

# 3. hand-pinned unit tests against the output file just built
SAF_ARTIFACT_DIR="${LOCK}" python -m pytest tests/unit/test_gse202379_saf_grades.py -v \
  > "${VAL}/unit_tests.txt" 2>&1

# 4. hold back: source.sha256, ARTIFACTS.json, COMPLETE
python scripts/recover_gse202379_saf_grades.py --finalize "${EXEC}"

echo "=== done $(date) ==="
echo "EXEC=${EXEC}"
