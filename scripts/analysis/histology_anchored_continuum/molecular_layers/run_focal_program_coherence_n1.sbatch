#!/usr/bin/env bash
#SBATCH --job-name=hac_coh_n1
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/program_coherence/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/83_focal_program_coherence_n1.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_COHERENCE_N1_OUT:?HAC_COHERENCE_N1_OUT is unset}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" "${HAC_COHERENCE_N1_OUT}"

sha256sum "${HAC_COHERENCE_N1_OUT}"/*.tsv
