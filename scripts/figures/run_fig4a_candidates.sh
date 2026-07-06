#!/bin/bash
#SBATCH --job-name=matplotlib
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4a_cand_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4a_cand_%j.err
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
export MPLCONFIGDIR=/scratch/claude-91825/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/9ae1b7d3-c694-44cb-b3d3-abdd560d9a18/scratchpad/mplcache_$SLURM_JOB_ID
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python scripts/figures/fig4a_overview_candidates.py
echo "FIG4A_CAND_JOB_DONE $(date)"
