#!/bin/bash
#SBATCH --job-name=S4_sccoda
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/S4_sccoda_%j.out
#SBATCH --error=logs/S4_sccoda_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

# Install scCODA if missing — use --no-deps to avoid rpy2 build failure on this
# cluster, then install only the required runtime deps (tensorflow-probability
# is pulled in automatically by the import path; rpy2 is only used by an
# optional R-interface inside sccoda that we don't invoke).
if ! python -c "import sccoda" 2>/dev/null; then
  echo "scCODA not installed; attempting --no-deps install"
  pip install --quiet --user --no-deps sccoda || true
  # Ensure runtime deps that come with rnaseq are sufficient; install absent
  # ones individually (no rpy2).
  for PKG in tensorflow-probability matplotlib seaborn arviz; do
    python -c "import ${PKG//-/_}" 2>/dev/null || pip install --quiet --user "$PKG" || true
  done
fi

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export ALLOW_SCCODA_MISSING=1

# Run with graceful skip-on-missing-scCODA semantics in the Python script
python Analysis/SingleCell/scripts/348_compositional_proportion_retest.py || {
  EC=$?
  if [[ $EC -eq 0 ]]; then exit 0; fi
  echo "[warn] 348 exited $EC; checking if it produced partial output"
  ls -la Analysis/SingleCell/results_gpu_v2/proportion_compositional/ 2>&1 || true
  # Don't fail the chain — propeller / S3 binomial GLM remain as compositional alternatives
  exit 0
}
