#!/bin/bash
#SBATCH --job-name=edgeR
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_fixedeffect_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_fixedeffect_%j.err

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

# STAR-native fixed-effect (cohort-as-fixed) sensitivity via edgeR quasi-likelihood.
# mega_validation arm 1 reads canonical -s 2 merged_dge.rds + meta_matched.rds,
# writes ql_results.csv + qc.csv (self-contained; no atlas touch).
echo "=== Fixed-effect cohort sensitivity (edgeR-QL, STAR -s 2) ==="
echo "Started: $(date)"
Rscript analysis/integration/scripts/mega_validation/01_edgeRql_cohortFE.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: edgeR-QL cohortFE"; exit 1; fi
echo "Finished: $(date)"
