#!/usr/bin/env bash
#SBATCH --job-name=borzoi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/borzoi_logsed_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/borzoi_logsed_%j.err

# ---------------------------------------------------------------------------
# DEFINITIVE on-recipe Borzoi logSED eQTL direction gate (Stage-2a seqfunc).
# 4-fold ensemble (replicate-0..3) x forward+RC x GTEx-liver tissue-matched
# logSED over the 234 Broadaway lead signals (+ SORT1 anchor), then the audited
# concordance benchmark. Additive-only (new files, script 64b). Uses the `borzoi`
# micromamba env for scoring and `spatial` for the benchmark.
# ---------------------------------------------------------------------------
set -uo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
SRC="${ROOT}/GWAS/finemapping/src"

FASTA=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
TARGETS="${SEQFUNC}/tracks/targets_human.txt"
LIVER_TRACKS="${SEQFUNC}/tracks/liver_tracks_selected.csv"
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
TRUTH="${SEQFUNC}/broadaway_benchmark_truth.tsv"
SUBSTRATE="${SEQFUNC}/variant_substrate_hg38.tsv"
OUT="${SEQFUNC}/borzoi_eqtl_logsed_scores.tsv"

mkdir -p "${SEQFUNC}/logs"

echo "=== Borzoi on-recipe logSED eQTL scorer (4-fold ensemble, fwd+RC) ==="
echo "Node: $(hostname) | Date: $(date)"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || true
echo ""

eval "$(micromamba shell hook -s bash)"

# --- Stage 1: GPU logSED scorer (borzoi env) --------------------------------
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

python3 "${SRC}/64b_borzoi_eqtl_logsed.py" \
    --fasta        "${FASTA}" \
    --targets      "${TARGETS}" \
    --liver-tracks "${LIVER_TRACKS}" \
    --gtf          "${GTF}" \
    --truth        "${TRUTH}" \
    --substrate    "${SUBSTRATE}" \
    --folds        johahi/borzoi-replicate-0 johahi/borzoi-replicate-1 \
                   johahi/borzoi-replicate-2 johahi/borzoi-replicate-3 \
    --out          "${OUT}"
SCORER_RC=$?
echo "scorer exit code: ${SCORER_RC} (nonzero = SORT1 anchor non-positive or no leads; benchmark still runs)"

# --- Stage 2: audited concordance benchmark (spatial env) -------------------
if [ ! -s "${OUT}" ]; then
    echo "FATAL: scores file ${OUT} missing/empty; cannot benchmark."
    exit 1
fi

micromamba deactivate || true
micromamba activate spatial
export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1
cd "${SRC}"
python3 "${SRC}/64b_concordance_benchmark_logsed.py"
BENCH_RC=$?

echo "Done: $(date) | scorer_rc=${SCORER_RC} bench_rc=${BENCH_RC}"
