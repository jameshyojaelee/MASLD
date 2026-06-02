#!/bin/bash
#SBATCH --job-name=NICHES_install
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=Analysis/Perturbation/scripts/d4_circuit/logs/install_niches_%j.out
#SBATCH --error=Analysis/Perturbation/scripts/d4_circuit/logs/install_niches_%j.err

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

ENV_PREFIX=.envs/perturbation_ccc

echo "[$(date)] Host: $(hostname)"
echo "[$(date)] Activating env at: $ENV_PREFIX"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate "$ENV_PREFIX"

echo "[$(date)] R version:"
R --version | head -1

RSCRIPT=/tmp/install_niches_${SLURM_JOB_ID:-$$}.R
cat > "$RSCRIPT" <<'RSCRIPT_EOF'
options(timeout = 900, repos = c(CRAN = "https://cloud.r-project.org"))

cat("--- Installing optparse + plyr from CRAN ---\n")
install.packages(c("optparse", "plyr"))

cat("--- Installing BisRNA (archived) from CRAN Archive ---\n")
if (!requireNamespace("BisRNA", quietly = TRUE)) {
  remotes::install_url(
    "https://cran.r-project.org/src/contrib/Archive/BisRNA/BisRNA_0.2.2.tar.gz",
    upgrade = "never"
  )
}
cat("BisRNA installed:", requireNamespace("BisRNA", quietly = TRUE), "\n")

cat("--- Installing NICHES from GitHub msraredon/NICHES@dev ---\n")
remotes::install_github("msraredon/NICHES", ref = "dev", upgrade = "never")

cat("--- Verifying ---\n")
suppressPackageStartupMessages({
  library(NICHES)
  library(optparse)
})
cat("NICHES version:", as.character(packageVersion("NICHES")), "\n")
cat("optparse version:", as.character(packageVersion("optparse")), "\n")
cat("OK\n")
RSCRIPT_EOF

echo "[$(date)] Running install R script: $RSCRIPT"
Rscript "$RSCRIPT"

echo "[$(date)] Verifying installs ..."
Rscript -e 'suppressPackageStartupMessages({library(NICHES); library(optparse)}); cat("NICHES version:", as.character(packageVersion("NICHES")), "\n"); cat("optparse version:", as.character(packageVersion("optparse")), "\n"); cat("OK\n")'

echo "[$(date)] DONE"
