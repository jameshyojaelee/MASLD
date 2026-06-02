#!/bin/bash
#SBATCH --job-name=liver_bulk_GSE167523
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=logs/bulk_%j.out
#SBATCH --error=logs/bulk_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# 1. Activate Environment
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
micromamba activate rnaseq

# 2. Define Directories
# Resolve Project Root (Robust to SLURM spooling)
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
PROJECT_ROOT=$(dirname "$SUBMIT_DIR")
cd "$PROJECT_ROOT"
echo "Running in PROJECT_ROOT: $PROJECT_ROOT"
export SNAKEMAKE_OUTPUT_CACHE="${PROJECT_ROOT}/.snakemake/cache"

echo "Starting BigMem Bulk Processing on $(hostname)"
echo "Allocated: 60 CPUs, 500GB RAM"

# 3. Unlock directory (just in case previous jobs left locks)
snakemake --unlock --cores 1 --configfile workflow/config.yaml -s workflow/Snakefile

# 4. Run Snakemake in Local Mode
# --cores 60: Uses the 60 allocated CPUs as local resources
# --resources: Inform Snakemake of limits (though local mode respects --cores)
# --keep-going: Continue other parallel branches if one sample fails
# --rerun-incomplete: Fix any files left by cancelled jobs
snakemake --cores 60 \
    --resources mem_mb=500000 \
    --configfile workflow/config.yaml \
    -s workflow/Snakefile \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --latency-wait 60

echo "Pipeline finished."
