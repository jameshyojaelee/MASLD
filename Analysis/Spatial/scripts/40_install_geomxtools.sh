#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --job-name=install_geomxtools
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/install_geomxtools_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/install_geomxtools_%j.err

# Create a dedicated GeoMx env via conda/bioconda (R 4.4 + GeomxTools 3.10.0).
# Source-installing GeomxTools into rnaseq fails because conda R's Makeconf
# has SHLIB_LIBADD empty, which breaks dotCall64 / roptim native code.
set -eo pipefail
set +u
eval "$(micromamba shell hook --shell bash)"

ENV_NAME=geomx
if micromamba env list | awk '{print $1}' | grep -q "^${ENV_NAME}$"; then
  echo "[$(date)] env ${ENV_NAME} already exists — updating"
  micromamba install -y -n ${ENV_NAME} \
    -c bioconda -c conda-forge \
    bioconductor-geomxtools=3.10.0 \
    bioconductor-nanostringnctools \
    r-readxl r-writexl r-dplyr r-tibble
else
  echo "[$(date)] Creating env ${ENV_NAME}"
  micromamba create -y -n ${ENV_NAME} \
    -c bioconda -c conda-forge \
    "r-base=4.4" \
    bioconductor-geomxtools=3.10.0 \
    bioconductor-nanostringnctools \
    r-readxl r-writexl r-dplyr r-tibble
fi

echo "[$(date)] Verifying installation"
micromamba activate ${ENV_NAME}
R -e 'library(NanoStringNCTools); library(GeomxTools); packageVersion("GeomxTools"); packageVersion("NanoStringNCTools"); sessionInfo()'
echo "[$(date)] DONE"
