#!/bin/bash
#SBATCH --job-name=power_saturation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --output=logs/power_saturation_%j.out
#SBATCH --error=logs/power_saturation_%j.err

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== Power saturation analysis ==="
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start: $(date)"
echo ""

Rscript analysis/integration/scripts/15c_power_saturation.R

echo ""
echo "End: $(date)"
