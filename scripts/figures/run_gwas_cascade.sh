#!/bin/bash
#SBATCH --job-name=matplotlib
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/cascade_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/cascade_%j.err
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
export MPLCONFIGDIR=/scratch/claude-91825/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/c00babb4-8b11-487d-8d31-b775efb6917b/scratchpad/mplcache_$SLURM_JOB_ID
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python scripts/figures/fig1_gwas_creative_options.py
echo "CASCADE_JOB_DONE $(date)"
