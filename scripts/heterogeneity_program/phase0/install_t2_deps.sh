#!/bin/bash
#SBATCH --job-name=deps
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --output=scripts/heterogeneity_program/phase0/logs/deps_%j.out
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/heterogeneity_program/phase0/logs
source ~/.bashrc
micromamba activate rnaseq

# Additive R-level installs (do not touch conda-managed core packages).
Rscript -e '
ok <- function(p) requireNamespace(p, quietly=TRUE)
cran <- c("mclust","archetypes","diptest","energy","infotheo","remotes")
for (p in cran) if (!ok(p)) try(install.packages(p, repos="https://cloud.r-project.org", quiet=TRUE))
if (!ok("BiocManager")) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
bioc <- c("miloR","speckle","ConsensusClusterPlus")
for (p in bioc) if (!ok(p)) try(BiocManager::install(p, update=FALSE, ask=FALSE))
try(remotes::install_github("vitkl/ParetoTI", upgrade="never", quiet=TRUE))
cat("\n── R dependency verification ──\n")
for (p in c("mclust","archetypes","diptest","energy","infotheo","miloR","speckle","ConsensusClusterPlus","ParetoTI"))
  cat(sprintf("  %-22s %s\n", p, if (ok(p)) "OK" else "MISSING"))
'
# scCODA (python side of rnaseq env)
python -c "import importlib.util as u; print('scCODA', 'OK' if u.find_spec('sccoda') else 'MISSING')" || true
pip install --quiet sccoda 2>&1 | tail -3 || echo "scCODA pip install attempted"
python -c "import importlib.util as u; print('scCODA after', 'OK' if u.find_spec('sccoda') else 'MISSING')" || true
echo "DONE deps"
