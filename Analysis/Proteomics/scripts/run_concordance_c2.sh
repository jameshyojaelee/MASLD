#!/bin/bash
#SBATCH --job-name=concordance
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Proteomics/logs/concordance_c2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Proteomics/logs/concordance_c2_%j.err

# Re-run proteomics concordance against C2 canonical bulk DEG (limma_voom_qw__C2).
# Both scripts already repointed to canonical_deg_results.csv (2026-06-08);
# this refreshes their downstream outputs (fgsea NES + stratified rho table).
set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PROTEO_DIR="${BASE}/Analysis/Proteomics"
mkdir -p "${PROTEO_DIR}/logs"

echo "=== Proteomics concordance vs C2 ==="
echo "Time: $(date)"; echo "Node: $(hostname)"

eval "$(micromamba shell hook -s bash)"
set +u
micromamba activate rnaseq
set -u

export MASLD_PROJECT_ROOT="${BASE}"

echo "--- [1/2] differential_proteomics.R (DE + concordance v3 + fgsea NES) ---"
Rscript "${PROTEO_DIR}/scripts/differential_proteomics.R"

echo "--- [2/2] mrna_protein_concordance.R (stratified rho table + figure) ---"
Rscript "${PROTEO_DIR}/scripts/mrna_protein_concordance.R"

echo "Done: $(date)"
