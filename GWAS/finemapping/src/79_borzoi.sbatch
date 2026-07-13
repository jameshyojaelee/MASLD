#!/usr/bin/env bash
#SBATCH --job-name=borzoi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_borzoi_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_borzoi_%j.err

# Phase-5 src/79 Borzoi caQTL ACCESSIBILITY corroboration: signed ref->alt delta over
# the liver ATAC (2035 hepatocyte / 2023) + DNASE (1302/1303/1367/1510/1623) tracks,
# over Currin caQTL LEADS. Modality-matched (accessibility tracks, NOT gene-expr).
# GPU-heavy (524kb x2 fwd passes) -> subsample via MAXVARS. ADDITIVE; edits nothing.
set -uo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC="${ROOT}/GWAS/finemapping/src"
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
mkdir -p "${SEQFUNC}/logs"

echo "=== src/79 Borzoi caQTL accessibility panel ==="; echo "$(hostname) | $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

eval "$(micromamba shell hook -s bash)"
micromamba activate borzoi
export PYTHONUNBUFFERED=1 PYTHONNOUSERSITE=1
export HF_HOME="${SEQFUNC}/hf_cache" HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
SITE_PKGS="${CONDA_PREFIX}/lib/python3.11/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

MAXVARS="${MAXVARS:-3000}"
python3 "${SRC}/79_score_panels.py" --panel borzoi --max-variants "${MAXVARS}"
echo "Done: $(date)"
