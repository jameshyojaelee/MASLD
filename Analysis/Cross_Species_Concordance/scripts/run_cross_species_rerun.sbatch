#!/bin/bash
#SBATCH --job-name=xspecies_rerun
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/cross_species_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/cross_species_rerun_%j.out
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts"

echo "=== Cross-Species Concordance Rerun (1,444 samples) ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK  Mem: $SLURM_MEM_PER_NODE"
echo ""

# 02b and 03a are independent: different inputs (both read merged_dge.rds but
# write to different output files) and can safely run in parallel.

echo "[02b] ssGSEA concordance starting in background..."
Rscript "$SCRIPT_DIR/02b_ssgsea_concordance.R" &
PID_02B=$!

echo "[03a] WGCNA module preservation starting in background..."
Rscript "$SCRIPT_DIR/03a_wgcna_module_preservation.R" &
PID_03A=$!

echo "Both scripts launched. Waiting for completion..."
echo ""

wait $PID_02B
RC_02B=$?
if [ $RC_02B -ne 0 ]; then
  echo "ERROR: 02b_ssgsea_concordance.R failed with exit code $RC_02B"
else
  echo "[02b] ssGSEA concordance COMPLETE"
fi

wait $PID_03A
RC_03A=$?
if [ $RC_03A -ne 0 ]; then
  echo "ERROR: 03a_wgcna_module_preservation.R failed with exit code $RC_03A"
else
  echo "[03a] WGCNA module preservation COMPLETE"
fi

echo ""
echo "=== Rerun Complete ==="
echo "Date: $(date)"
echo ""
echo "Results:"
ls -lh /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/results/

if [ $RC_02B -ne 0 ] || [ $RC_03A -ne 0 ]; then
  exit 1
fi
