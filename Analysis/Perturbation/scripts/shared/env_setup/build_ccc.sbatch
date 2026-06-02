#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=8:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=env_ccc
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_ccc.log

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Prevent user-site (~/.local/lib/pythonX.Y/site-packages) leak into env Python.
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

ENV_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_ccc
YAML=envs/perturbation/perturbation_ccc.yml
export LD_LIBRARY_PATH="$ENV_PATH/lib:${LD_LIBRARY_PATH:-}"

echo "=========================================="
echo "Building $ENV_PATH from $YAML"
echo "Start time: $(date -Iseconds)"
echo "NOTE: NICHES is R-only; installed via R remotes post-solve"
echo "=========================================="

if [ -d "$ENV_PATH" ]; then
  echo "Removing existing partial env at $ENV_PATH"
  rm -rf "$ENV_PATH"
fi

# Patch YAML to strip cellchat-py placeholder (no PyPI package guaranteed)
PATCHED_YAML=envs/perturbation/perturbation_ccc_patched.yml
python3 - <<PY
import re
src = open("$YAML").read()
src = re.sub(r"^\s*-\s*cellchat-py.*$", "    # cellchat-py removed (no canonical PyPI pkg; CellChat is R)", src, flags=re.M)
open("$PATCHED_YAML","w").write(src)
print("Patched YAML written:", "$PATCHED_YAML")
PY

T0=$(date +%s)
micromamba env create -f "$PATCHED_YAML" -p "$ENV_PATH" -y || {
  echo "ERROR: micromamba env create failed for CCC env"
  exit 1
}
T1=$(date +%s)
echo "Solve+install time: $((T1-T0))s"

echo "=== Install NICHES via R remotes ==="
micromamba run -p "$ENV_PATH" R --no-save <<'RSCRIPT'
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
tryCatch({
  remotes::install_github("msraredon/NICHES", ref = "master", upgrade = "never")
  cat("NICHES install ATTEMPT 1 (ref=master) result: ", as.character(requireNamespace("NICHES", quietly = TRUE)), "\n")
}, error = function(e) {
  cat("ATTEMPT 1 failed:", conditionMessage(e), "\n")
  tryCatch({
    remotes::install_github("msraredon/NICHES", ref = "dev", upgrade = "never")
    cat("NICHES install ATTEMPT 2 (ref=dev) success:", as.character(requireNamespace("NICHES", quietly = TRUE)), "\n")
  }, error = function(e2) cat("ATTEMPT 2 also failed:", conditionMessage(e2), "\n"))
})
sessionInfo()
RSCRIPT

echo "=== Verification ==="
micromamba run -p "$ENV_PATH" python -c "
import sys, importlib
print('Python:', sys.version)
mods = ['torch','scanpy','anndata','squidpy','omnipath','ot','liana','commot','pertpy']
for m in mods:
    try:
        mod = importlib.import_module(m)
        v = getattr(mod, '__version__', '?')
        print(f'  OK {m} {v}')
    except Exception as e:
        print(f'  FAIL {m}: {e}')
try:
    import torch
    print('  CUDA available:', torch.cuda.is_available())
except Exception:
    pass
"

micromamba run -p "$ENV_PATH" R --no-save -e '
cat("R version:", R.version$version.string, "\n")
for (p in c("Seurat","NICHES","tidyverse","SingleCellExperiment")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %s %s: %s\n", if (ok) "OK" else "FAIL", p,
              if (ok) as.character(packageVersion(p)) else "missing"))
}
'
echo "End time: $(date -Iseconds)"
echo "DONE: CCC env build"
