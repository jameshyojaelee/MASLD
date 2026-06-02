#!/bin/bash
#SBATCH --job-name=liver_progression
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/progression_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/progression_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlee@nygenome.org

set -euo pipefail

# Activate R environment
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "============================================================"
echo "  CONTINUOUS DISEASE PROGRESSION MODELING"
echo "  Date: $(date)"
echo "  Job ID: $SLURM_JOB_ID"
echo "============================================================"

# Script 13: NAFL vs NASH
echo -e "\n=== 13: NAFL vs NASH DE ==="
Rscript scripts/13_nafl_vs_nash_de.R

# Script 14: Fibrosis Progression
echo -e "\n=== 14: Fibrosis Progression DE ==="
Rscript scripts/14_fibrosis_progression_de.R

# Script 15: NAS Component DE
echo -e "\n=== 15: NAS Component DE ==="
Rscript scripts/15_nas_component_de.R

# Script 16: Consensus Signatures
echo -e "\n=== 16: Consensus Signatures ==="
Rscript scripts/16_disease_signatures_consensus.R

echo "============================================================"
echo "  PROGRESSION MODELING COMPLETE"
echo "  Finished: $(date)"
echo "============================================================"
