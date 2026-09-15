#!/usr/bin/env bash
#SBATCH --job-name=agp_c1-endpoint2-v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2v2_%j.err
set -uo pipefail

# C1 endpoint 2 as the FOUR-model comparison the spec describes, now that the two Borzoi GPU jobs
# (21764806, 21764807) have landed with 3,102 variants scored and 0 unresolved. No GPU work here;
# the scores exist. Reads the GPU outputs from the c1-endpoint2-gpu package read-only and writes
# only into this package.

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${ROOT}/GWAS/finemapping/results/alphagenome_program/c1-endpoint2-v2-20260915T084559Z
SCORED=${ROOT}/GWAS/finemapping/results/alphagenome_program/c1-endpoint2-gpu-20260914T203750Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

echo "host $(hostname) start $(date -u)"
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8

"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/c1e2_05_recompute_v2.py" \
  --outdir "${OUT}" --scoredir "${SCORED}"
rc=$?
echo "recompute exit ${rc}"
[ ${rc} -ne 0 ] && exit ${rc}

"${PY}" -m pip freeze > "${OUT}/pip_freeze.txt" 2>/dev/null
"${PY}" -c "import sys; print(sys.version)" > "${OUT}/python_version.txt"
cp "${ROOT}/scripts/analysis/alphagenome_program/c1e2_05_recompute_v2.py" \
   "${ROOT}/scripts/analysis/alphagenome_program/c1e2_05_recompute_v2.sbatch" \
   "${OUT}/sbatch/"
( cd "${OUT}/sbatch" && sha256sum ./* > "${OUT}/analysis_scripts.sha256" )
cat > "${OUT}/seeds.json" <<JSON
{
  "block_bootstrap_seed": 20260914,
  "block_bootstrap_draws": 10000,
  "resampling_unit": "1-Mb block, chrN~Mb",
  "borzoi_model": "johahi/borzoi-replicate-0",
  "borzoi_seq_len_bp": 524288,
  "borzoi_window_bins_either_side": 16,
  "borzoi_dtype": "float32",
  "chrombpnet_scorer_seed": 42,
  "chrombpnet_batch_size": 1024,
  "note": "no seed in this package differs from c1-endpoint2-gpu-20260914T203750Z/seeds.json"
}
JSON
echo "done $(date -u)"
