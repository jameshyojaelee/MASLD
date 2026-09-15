#!/usr/bin/env bash
#SBATCH --job-name=pxd_two_axis
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/slurm-%x-%j.out

set -euo pipefail
BENCH="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
SCRIPT="${BENCH}/scripts/run_pxd051911_two_axis_protein.py"
PY="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python3"
: "${PXD_TWO_AXIS_OUT:?PXD_TWO_AXIS_OUT is unset}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${PY}" -c "import ast,sys; ast.parse(open('${SCRIPT}').read())"
"${PY}" "${SCRIPT}" "${PXD_TWO_AXIS_OUT}"

sha256sum "${PXD_TWO_AXIS_OUT}"/*.tsv "${PXD_TWO_AXIS_OUT}"/*.txt
