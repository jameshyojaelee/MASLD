#!/bin/bash
#SBATCH --job-name=strandcert
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/strand_cert/cert_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/strand_cert/cert_%j.err

# A real script rather than --wrap: sbatch --wrap executes under /bin/sh here,
# and `eval "$(micromamba shell hook -s bash)"` is bash-only syntax.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript /scratch/claude-91825/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/493bb570-0e68-4a71-813d-0caa8dadd59b/scratchpad/test_certificates.R
