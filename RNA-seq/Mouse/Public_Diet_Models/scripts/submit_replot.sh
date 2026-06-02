#!/bin/bash
#SBATCH --job-name=replot_volcanos
#SBATCH --output=logs/replot_volcanos_%j.out
#SBATCH --error=logs/replot_volcanos_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq

MAMBA_ENV="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/unified_env"

echo "Starting volcano replotting at $(date)"
micromamba run -p $MAMBA_ENV Rscript scripts/replot_volcanos_publication.R
echo "Completed at $(date)"
