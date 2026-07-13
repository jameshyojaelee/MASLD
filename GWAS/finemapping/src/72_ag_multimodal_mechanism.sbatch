#!/usr/bin/env bash
#SBATCH --job-name=alphagenome
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/alphagenome_mech_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/alphagenome_mech_%j.err

# ---------------------------------------------------------------------------
# Phase-4 seqfunc: AlphaGenome MULTIMODAL MECHANISM profile (script 72).
# Full head set (contact map + histone/TF ChIP + CAGE/PRO-cap + polyA + splice +
# accessibility) -> per-variant mechanism-CLASS / which-gene / magnitude profile.
# Additive-only. NO GPU: AlphaGenome inference is Google-side; this job only opens
# a gRPC channel to the AlphaGenome endpoint. Runs in the `alphagenome` env
# (python 3.11 + alphagenome 0.7.0).
#
# Targets: 81 eQTL-absent regulatory leads + 42 coding effectors (30 in-target) +
# SORT1/HSD17B13/PNPLA3 anchors.  API key: $ALPHA_GENOME_API_KEY or ~/.alphagenome_key
# (chmod 600); NEVER logged or committed.  Args passed through (e.g. --anchors-only,
# --limit N, --dry-run) for smoke tests.
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC="${ROOT}/GWAS/finemapping/src"
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
mkdir -p "${SEQFUNC}/logs"

echo "=== AlphaGenome multimodal mechanism (script 72) ==="
echo "Node: $(hostname) | Date: $(date)"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate alphagenome
export PYTHONUNBUFFERED=1
export PYTHONNOUSERSITE=1

python3 "${SRC}/72_ag_multimodal_mechanism.py" \
    --out "${SEQFUNC}/ag_mechanism_profile.tsv" \
    "$@"

echo "Done: $(date)"
