#!/bin/bash
#SBATCH --job-name=liver_megabulk
#SBATCH --output=logs/megabulk_%j.out
#SBATCH --error=logs/megabulk_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=70
#SBATCH --mem=1500GB
#SBATCH --time=48:00:00

echo "ERROR: PRJNA512027 is retired/noncanonical; pipeline submission is prohibited." >&2
exit 64

set -euo pipefail

# 1. Activate Environment
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq
set -u

# 2. Navigate to Pipeline Directory
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027

# 3. Unlock directory (good practice for restarts)
snakemake --unlock --cores 1

# 4. Run Snakemake LOCALLY on this BigMem node
# We use --cores 70 to utilize the full allocation of this job.
# No --profile needed because we are NOT submitting jobs to SLURM; we are the job.
echo "Starting Mega Bulk Processing for PRJNA512027 on $(hostname)"
echo "Allocated: 70 CPUs, 1500GB RAM"

snakemake \
    --cores 70 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --latency-wait 60 \
    --configfile workflow/config.yaml
