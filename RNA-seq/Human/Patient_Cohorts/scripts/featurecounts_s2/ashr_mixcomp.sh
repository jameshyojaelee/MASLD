#!/bin/bash
#SBATCH --job-name=ashr_mixcomp
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/logs/ashr_mixcomp_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/logs/ashr_mixcomp_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== ashr mixcomp sensitivity (halfuniform vs normal, -s 2). Writes audit_sensitivity/ashr_mixcomp_sensitivity/ ==="
echo "Started: $(date)  Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/05b2_ashr_normal_sensitivity.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 05b2_ashr_normal_sensitivity.R"; exit 1; fi

echo "Finished: $(date)"
