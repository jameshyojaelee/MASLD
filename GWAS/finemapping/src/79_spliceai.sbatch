#!/usr/bin/env bash
#SBATCH --job-name=spliceai
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_spliceai_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_spliceai_%j.err

# Phase-5 src/79 SpliceAI panel: DS_AG/AL/DG/DL + max_ds over sQTL LEADS
# (leakage_exploratory, magnitude-only). CPU (CUDA disabled). ADDITIVE; edits nothing.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC="${ROOT}/GWAS/finemapping/src"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
set +u; micromamba activate "${ROOT}/.envs/spliceai_throwaway"; set -u
export CUDA_VISIBLE_DEVICES="" TF_CPP_MIN_LOG_LEVEL=3 PYTHONUNBUFFERED=1
MAXVARS="${MAXVARS:-6000}"   # sQTL exploratory: cap the 137k lead-introns (top by p)
python "${SRC}/79_score_panels.py" --panel spliceai --max-variants "${MAXVARS}"
echo "Done: $(date)"
