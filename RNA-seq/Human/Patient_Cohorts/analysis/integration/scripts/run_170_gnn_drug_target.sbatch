#!/bin/bash
#SBATCH --job-name=gnn_drug_target
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/170_gnn_drug_target_%j.out
#SBATCH --error=logs/170_gnn_drug_target_%j.err

set -euo pipefail

echo "=== Script 170: GNN Drug-Target Predictor ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

# Use spatial env (has torch, torch_geometric, sklearn without ABI conflicts)
eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

# Set thread counts
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p logs results/novel_ml/gnn_drug_target

python scripts/170_gnn_drug_target.py

echo ""
echo "Finish: $(date)"
echo "=== Done ==="
