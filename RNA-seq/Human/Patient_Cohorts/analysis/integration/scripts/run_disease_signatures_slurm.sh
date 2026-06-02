#!/bin/bash
#SBATCH --job-name=liver_disease_sigs
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/disease_sigs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/disease_sigs_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=cpu

set -euo pipefail

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"

echo "=== Disease Signature Analyses: $(date) ==="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK"

# Activate environment
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Script 13: NAFL vs NASH
echo ""
echo "=== 13: NAFL vs NASH DE ==="
Rscript "${SCRIPTS}/13_nafl_vs_nash_de.R"

# Script 14: Fibrosis progression
echo ""
echo "=== 14: Fibrosis Progression DE ==="
Rscript "${SCRIPTS}/14_fibrosis_progression_de.R"

# Script 15: NAS components
echo ""
echo "=== 15: NAS Component DE ==="
Rscript "${SCRIPTS}/15_nas_component_de.R"

# Script 16: Consensus integration
echo ""
echo "=== 16: Disease Signatures Consensus ==="
Rscript "${SCRIPTS}/16_disease_signatures_consensus.R"

echo ""
echo "=== All disease signature scripts completed: $(date) ==="
