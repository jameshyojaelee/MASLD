#!/usr/bin/env bash
#SBATCH --job-name=agp_c1-endpoint2-gpu-bzrepro
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=b6k
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_bzrepro_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_bzrepro_%j.err
set -uo pipefail

# Rescore 200 Currin leads that the archived src/79 Borzoi panel already scored, with this
# package's loader, so the mixed-provenance Borzoi column carries a reproduction check.

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PKG=${ROOT}/GWAS/finemapping/results/alphagenome_program/c1-endpoint2-gpu-20260914T203750Z
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"

echo "host $(hostname) start $(date -u)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

export PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "${ROOT}/scripts/analysis/alphagenome_program/c1e2_04_build_repro_list.py" \
  --out "${PKG}/scoring/borzoi/borzoi_repro_targets.tsv" --n 200

eval "$(micromamba shell hook -s bash)"
micromamba activate borzoi
export PYTHONUNBUFFERED=1
export HF_HOME="${SEQFUNC}/hf_cache" HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
SITE_PKGS="${CONDA_PREFIX}/lib/python3.11/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

python3 "${ROOT}/scripts/analysis/alphagenome_program/c1e2_02_borzoi_score.py" \
  --targets "${PKG}/scoring/borzoi/borzoi_repro_targets.tsv" \
  --out "${PKG}/scoring/borzoi/borzoi_repro_scores.tsv" \
  --checkpoint-every 50
echo "done $(date -u)"
