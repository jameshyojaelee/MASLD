#!/usr/bin/env bash
#SBATCH --job-name=borzoi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/borzoi_magnitude_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/borzoi_magnitude_%j.err

# ---------------------------------------------------------------------------
# AXIS 4: Borzoi MAGNITUDE nomination + gene attribution (direction-INDEPENDENT).
# For each eQTL-ABSENT fine-mapped signal lead, enumerate protein_coding/lncRNA
# candidate genes in the Borzoi window and score each with the 64b logSED recipe
# (4-fold ensemble x fwd+RC x GTEx-liver tissue-matched); report best-magnitude
# gene + magnitude percentile vs the eQTL-present lead distribution. ADDITIVE-only
# (new scripts 69a/69); imports 64b machinery verbatim, edits nothing.
# ---------------------------------------------------------------------------
set -uo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
SRC="${ROOT}/GWAS/finemapping/src"

FASTA=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
LIVER_TRACKS="${SEQFUNC}/tracks/liver_tracks_selected.csv"
LEADS="${SEQFUNC}/borzoi_magnitude_leads.tsv"
REFERENCE="${SEQFUNC}/borzoi_eqtl_logsed_scores.tsv"
OUT="${SEQFUNC}/borzoi_magnitude_nominations.tsv"
OUT_CAND="${SEQFUNC}/borzoi_magnitude_candidates.tsv"

mkdir -p "${SEQFUNC}/logs"

echo "=== Borzoi magnitude nomination (Axis 4, 4-fold ensemble, fwd+RC) ==="
echo "Node: $(hostname) | Date: $(date)"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || true
echo ""

# conda/micromamba activation scripts reference unbound vars (ADDR2LINE, ...),
# which trip `set -u`; guard every activate/deactivate with `set +u`.
set +u
eval "$(micromamba shell hook -s bash)"
set -u

# --- Stage 0: CPU lead selection (rnaseq env; deterministic, cheap) ---------
set +u; micromamba activate rnaseq; set -u
python3 "${SRC}/69a_borzoi_magnitude_prep.py" || { echo "FATAL: 69a prep failed"; exit 1; }
set +u; micromamba deactivate || true; set -u

# --- Stage 1: GPU candidate scoring (borzoi env) ----------------------------
set +u; micromamba activate borzoi; set -u
export PYTHONUNBUFFERED=1
export PYTHONNOUSERSITE=1
export HF_HOME="${SEQFUNC}/hf_cache"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1

SITE_PKGS="${CONDA_PREFIX}/lib/python3.11/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

python3 "${SRC}/69_borzoi_magnitude_nom.py" \
    --fasta          "${FASTA}" \
    --gtf            "${GTF}" \
    --liver-tracks   "${LIVER_TRACKS}" \
    --leads          "${LEADS}" \
    --reference      "${REFERENCE}" \
    --folds          johahi/borzoi-replicate-0 johahi/borzoi-replicate-1 \
                     johahi/borzoi-replicate-2 johahi/borzoi-replicate-3 \
    --out            "${OUT}" \
    --out-candidates "${OUT_CAND}"
RC=$?

echo "Done: $(date) | scorer_rc=${RC}"
exit ${RC}
