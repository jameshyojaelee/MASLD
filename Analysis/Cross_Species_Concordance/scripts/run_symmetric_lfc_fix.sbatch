#!/bin/bash
#SBATCH --job-name=xspecies_lfc_fix
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts/logs/symmetric_lfc_fix_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts/logs/symmetric_lfc_fix_%j.out
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts"

echo "=== Symmetric logFC Fix: Cross-Species Concordance ==="
echo "Fix: use raw dream logFC (dream_results.csv) instead of ashr shrunk logFC"
echo "      (dream_results_ashr.csv) for human disease-vs-control signature."
echo "      Mouse side already uses raw logFC. This eliminates asymmetric shrinkage bias."
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK  Mem: $SLURM_MEM_PER_NODE"
echo ""

echo "[01] Gene-level concordance (4x5 matrix)..."
Rscript "$SCRIPT_DIR/01_corrected_gene_concordance.R"
echo "[01] COMPLETE"
echo ""

echo "=== Fix Complete ==="
echo "Date: $(date)"
echo ""
echo "Updated results:"
ls -lh /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/results/gene_concordance_*.csv
