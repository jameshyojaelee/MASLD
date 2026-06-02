#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --output=logs/fm_%x_%A_%a.out
#SBATCH --error=logs/fm_%x_%A_%a.err

# Fine-mapping per-locus SLURM wrapper
# Called as array job: --array=1-N where N = number of lead SNPs
# Usage: sbatch --array=1-N 03_run_fm_per_locus.sh <sumstats_name> <ld_pop> <lead_snps_file> <N_tot> <N_cases> <window_mb> <ancestry>

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping

# Activate rnaseq environment
eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping

sumstats_name=$1
ld_pop=$2
lead_snps_file=$3
N_tot=$4
N_cases=$5
window_mb=$6
ancestry=$7

# Get the locus for this array task (skip header line)
LOCUS=$(awk -v line=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==line {print $3}' ${lead_snps_file})

if [ -z "$LOCUS" ]; then
    echo "ERROR: No locus found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Array task ${SLURM_ARRAY_TASK_ID}: Running fine-mapping for locus ${LOCUS}"
echo "Study: ${sumstats_name}, LD: ${ld_pop}, Ancestry: ${ancestry}"

Rscript src/03_run_fm_per_locus.R \
    ${sumstats_name} \
    ${ld_pop} \
    ${LOCUS} \
    ${N_tot} \
    ${N_cases} \
    ${window_mb} \
    ${ancestry}

echo "Completed locus ${LOCUS} (exit code: $?)"
