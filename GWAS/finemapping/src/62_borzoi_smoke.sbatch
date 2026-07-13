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
# Stage-2 seqfunc GPU smoke test: prove Borzoi runs on this cluster's GPU and
# scores a variant ref-vs-alt delta on liver RNA + hepatocyte/HepG2 accessibility
# tracks. Additive-only (new files, script number 62). Uses the throwaway
# `borzoi` micromamba env (borzoi-pytorch 0.5.1 + torch 2.13.0+cu130).
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
SRC="${ROOT}/GWAS/finemapping/src"

FASTA=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
TARGETS="${SEQFUNC}/tracks/targets_human.txt"
LIVER_TRACKS="${SEQFUNC}/tracks/liver_tracks_selected.csv"

echo "=== Borzoi GPU smoke test ==="
echo "Node: $(hostname) | Date: $(date)"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || true
echo ""

# --- env (isolated from ~/.local user-site, which shadows torch deps) --------
eval "$(micromamba shell hook -s bash)"
micromamba activate borzoi

export PYTHONUNBUFFERED=1
export PYTHONNOUSERSITE=1                      # ignore ~/.local (leaks torch 2.1 etc.)
export HF_HOME="${SEQFUNC}/hf_cache"           # shared gpfs cache (populated on login node)
export HF_HUB_OFFLINE=1                         # GPU nodes may be air-gapped -> use cache only
export TRANSFORMERS_OFFLINE=1

# bundled NVIDIA cu13 libs (belt-and-suspenders; torch usually resolves via RPATH)
SITE_PKGS="${CONDA_PREFIX}/lib/python3.11/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

python3 "${SRC}/62_borzoi_score.py" \
    --fasta        "${FASTA}" \
    --targets      "${TARGETS}" \
    --liver-tracks "${LIVER_TRACKS}" \
    --model        "johahi/borzoi-replicate-0" \
    --out-dir      "${SEQFUNC}" \
    --win-bins     16

echo "Done: $(date)"
