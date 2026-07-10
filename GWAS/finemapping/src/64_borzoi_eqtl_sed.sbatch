#!/usr/bin/env bash
#SBATCH --job-name=borzoi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/borzoi_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/borzoi_%j.err

# ---------------------------------------------------------------------------
# Stage-2a seqfunc: Borzoi GENE-LEVEL SED eQTL scorer over the 234 independent
# Broadaway lead signals (+ SORT1 positive-control anchor). Additive-only
# (new files, script number 64). Uses the `borzoi` micromamba env
# (borzoi-pytorch 0.5.1 + torch 2.13.0+cu130).
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
SRC="${ROOT}/GWAS/finemapping/src"

FASTA=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
TARGETS="${SEQFUNC}/tracks/targets_human.txt"
LIVER_TRACKS="${SEQFUNC}/tracks/liver_tracks_selected.csv"
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
TRUTH="${SEQFUNC}/broadaway_benchmark_truth.tsv"
SUBSTRATE="${SEQFUNC}/variant_substrate_hg38.tsv"
OUT="${SEQFUNC}/borzoi_eqtl_scores.tsv"

mkdir -p "${SEQFUNC}/logs"

echo "=== Borzoi gene-level SED eQTL scorer ==="
echo "Node: $(hostname) | Date: $(date)"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || true
echo ""

# --- env (isolated from ~/.local user-site) ---------------------------------
eval "$(micromamba shell hook -s bash)"
micromamba activate borzoi

export PYTHONUNBUFFERED=1
export PYTHONNOUSERSITE=1
export HF_HOME="${SEQFUNC}/hf_cache"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1

SITE_PKGS="${CONDA_PREFIX}/lib/python3.11/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

python3 "${SRC}/64_borzoi_eqtl_sed.py" \
    --fasta        "${FASTA}" \
    --targets      "${TARGETS}" \
    --liver-tracks "${LIVER_TRACKS}" \
    --gtf          "${GTF}" \
    --truth        "${TRUTH}" \
    --substrate    "${SUBSTRATE}" \
    --model        "johahi/borzoi-replicate-0" \
    --out          "${OUT}"

echo "Done: $(date)"
