#!/bin/bash
#SBATCH --job-name=mcp_loo_244832
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --mem=400G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/470b_loo244832_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/470b_loo244832_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"
source "$HOME/.bashrc"
micromamba activate rapids_singlecell

echo "[LOO] building LOO inputs for top-3 datasets..."
python 470b_build_loo_GSE244832.py

for DS in GSE244832 GSE202379 GSE185477; do
  NAME="global_loo_${DS}"
  echo "=== LOO run: ${NAME} ==="
  python 420_run_cnmf.py prepare \
    --input "atlas_cnmf_global_loo_${DS}.h5ad" --name "${NAME}" \
    --k "16" --n-iter 100 --num-hvg 2500
  export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
  python 420_run_cnmf.py factorize --name "${NAME}" --worker-index 0 --total-workers 1
  python 420_run_cnmf.py consensus --name "${NAME}" --k "16" --density-threshold 0.03 || true
done
echo "[LOO] DONE."
