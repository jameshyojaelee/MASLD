#!/bin/bash
#SBATCH --job-name=liver_integration
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/full_5cohort_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/full_5cohort_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=6:00:00

# Full 5-cohort integration pipeline (03-08) with GSE126848 included.
# Scripts 01+02 were already re-run interactively with GSE126848 (608 QC-passing samples).
# This job re-runs the integration, variance partition, dream,
# consensus, and pathway scripts on the updated data.

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== Job started: $(date) ==="
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
echo "Partition: bigmem | RAM: 400G | CPUs: 32"
echo ""

# Script 03: Integrate counts (re-merge with 608 QC-passing samples)
echo "=== 03: Integrate Counts ==="
Rscript analysis/integration/scripts/03_integrate_counts.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 03_integrate_counts.R"; exit 1; fi

# Script 04: Variance partition
echo ""
echo "=== 04: Variance Partition ==="
Rscript analysis/integration/scripts/04_variance_partition.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 04_variance_partition.R"; exit 1; fi

# Script 05: Dream mega-analysis (biggest bottleneck — uses $SLURM_CPUS_PER_TASK=32)
echo ""
echo "=== 05: Dream Mega-Analysis (32 CPUs) ==="
Rscript analysis/integration/scripts/05_dream_mega_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 05_dream_mega_analysis.R"; exit 1; fi

# Script 07: Consensus DEGs
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
echo "=== All 5-cohort integration scripts completed: $(date) ==="
