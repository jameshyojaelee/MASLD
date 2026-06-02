#!/bin/bash
#SBATCH --job-name=liver_integration_main
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/integration_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/integration_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=4:00:00

# Activate environment
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== Job started: $(date) ==="
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo ""

# Script 05: dream mega-analysis (uses SLURM_CPUS_PER_TASK for parallelization)
echo "=== 05: Dream Mega-Analysis ==="
Rscript analysis/integration/scripts/05_dream_mega_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 05_dream_mega_analysis.R"; exit 1; fi

# Script 07: Consensus DEGs (depends on dream + meta results)
echo ""
echo "=== 07: Consensus DEGs ==="
Rscript analysis/integration/scripts/07_consensus_degs.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 07_consensus_degs.R"; exit 1; fi

# Script 08: Pathway analysis
echo ""
echo "=== 08: Pathway Analysis ==="
Rscript analysis/integration/scripts/08_pathway_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 08_pathway_analysis.R"; exit 1; fi

echo ""
echo "=== All scripts completed: $(date) ==="
