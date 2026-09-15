#!/usr/bin/env bash
#SBATCH --job-name=agp_c1-endpoint2-v2-qc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2v2qc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2v2qc_%j.err
set -uo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${ROOT}/GWAS/finemapping/results/alphagenome_program/c1-endpoint2-v2-20260915T084559Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
echo "host $(hostname) start $(date -u)"
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/c1e2_06_qc_addendum.py" --outdir "${OUT}"
rc=$?
echo "qc exit ${rc}"
[ ${rc} -ne 0 ] && exit ${rc}
cp "${ROOT}/scripts/analysis/alphagenome_program/c1e2_06_qc_addendum.py" \
   "${ROOT}/scripts/analysis/alphagenome_program/c1e2_06_qc_addendum.sbatch" "${OUT}/sbatch/"
( cd "${OUT}/sbatch" && sha256sum ./* > "${OUT}/analysis_scripts.sha256" )
echo "done $(date -u)"
