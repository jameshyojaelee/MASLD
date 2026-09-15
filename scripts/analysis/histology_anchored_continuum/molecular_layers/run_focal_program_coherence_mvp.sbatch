#!/usr/bin/env bash
#SBATCH --job-name=hac_coh_mvp
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/program_coherence/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/81_focal_program_coherence_mvp.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_COHERENCE_MVP_OUT:?HAC_COHERENCE_MVP_OUT is unset}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" "${HAC_COHERENCE_MVP_OUT}"

sha256sum "${HAC_COHERENCE_MVP_OUT}"/*.tsv
