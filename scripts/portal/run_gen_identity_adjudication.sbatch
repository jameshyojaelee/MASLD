#!/bin/bash
#SBATCH --job-name=portal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/gen_identity_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/gen_identity_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PRODUCER="$PROJECT_ROOT/scripts/portal/adjudicate_gen_ensembl_identity.py"
SOURCE_DEFAULT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context/frozen_evidence_classes.tsv"
IDENTITY_DEFAULT="$PROJECT_ROOT/data/gencode_v49_gene_metadata.tsv.gz"
OUTPUT_DEFAULT="$PROJECT_ROOT/RNA-seq/results/evidence_passports/candidates/program-context-v2-candidate-2026-08-07/preflight/gen_identity_adjudication_v1"

PASSPORT_GEN_SOURCE="${PASSPORT_GEN_SOURCE:-$SOURCE_DEFAULT}"
PASSPORT_GENCODE_IDENTITY="${PASSPORT_GENCODE_IDENTITY:-$IDENTITY_DEFAULT}"
PASSPORT_GEN_IDENTITY_OUTPUT="${PASSPORT_GEN_IDENTITY_OUTPUT:-$OUTPUT_DEFAULT}"
PASSPORT_ANALYSIS_RELEASE_ID="${PASSPORT_ANALYSIS_RELEASE_ID:-program-context-v2-candidate-2026-08-07-passports-v1}"

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

micromamba run -n spatial python -B "$PRODUCER" \
  --source "$PASSPORT_GEN_SOURCE" \
  --identity "$PASSPORT_GENCODE_IDENTITY" \
  --output-root "$PASSPORT_GEN_IDENTITY_OUTPUT" \
  --release-id "$PASSPORT_ANALYSIS_RELEASE_ID"
