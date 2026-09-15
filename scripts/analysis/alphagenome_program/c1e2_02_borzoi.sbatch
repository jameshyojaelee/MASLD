#!/usr/bin/env bash
#SBATCH --job-name=agp_c1-endpoint2-gpu-borzoi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=b6k
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=01:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_borzoi_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_borzoi_%j.err
set -uo pipefail

# Borzoi (johahi port) liver ATAC and DNase allele scores for the 500 union allelic-imbalance
# sites and the 2,602 Currin peak leads that carry an archived AlphaGenome score but no Borzoi
# score (priority 0 and 1 of borzoi_targets.tsv). Those 3,102 are exactly what the four-model
# matched sets need; the archived AlphaGenome panel caps endpoint 1 at 3,845 leads, so the
# remaining 28,458 unscored leads cannot enlarge any matched set. The archived src/79 panel ran
# 3,000 variants in 16 minutes wall, so 1.5 h is ~5x buffer and backfills where a 12 h ask did
# not. The b6k constraint the assignment names is dropped: every b6k GPU (2 per node on 11 nodes)
# is allocated and both a 12 h ask (job 21764740) and a 1.5 h ask (job 21764753) were given the
# same start estimate two days out, while four L40S nodes each had 2 of 3 GPUs free. The scorer
# now runs float32 only and raises on out-of-memory rather than falling back to bfloat16, so a
# smaller GPU cannot silently put part of the Borzoi column on a different numeric path from the
# archived part. Same recipe as the archived panel otherwise; writes only into the package.

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PKG=${ROOT}/GWAS/finemapping/results/alphagenome_program/c1-endpoint2-gpu-20260914T203750Z
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"

echo "host $(hostname) start $(date -u)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

eval "$(micromamba shell hook -s bash)"
micromamba activate borzoi
export PYTHONUNBUFFERED=1 PYTHONNOUSERSITE=1
export HF_HOME="${SEQFUNC}/hf_cache" HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
SITE_PKGS="${CONDA_PREFIX}/lib/python3.11/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

python3 "${ROOT}/scripts/analysis/alphagenome_program/c1e2_02_borzoi_score.py" \
  --targets "${PKG}/scoring/borzoi/borzoi_targets.tsv" \
  --out "${PKG}/scoring/borzoi/borzoi_scores.tsv" \
  --max-variants 3102 \
  --checkpoint-every 250
echo "done $(date -u)"
