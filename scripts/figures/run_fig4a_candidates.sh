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

# MAIN render: Tier-1/2 (liver-specific) COLOC universe -> fig4a_overview_*.pdf
echo "=== MAIN: Tier-1/2 COLOC (placement==main; coloc 862 / universe 13,390) ==="
python scripts/figures/fig4a_overview_candidates.py

# SUPP/sensitivity render: full 50-GWAS COLOC universe (adds Tier-3/4). Distinct
# _supp_full50gwas filenames, so it NEVER overwrites the main PDFs.
echo "=== SUPP: full 50-GWAS COLOC (coloc 1,234 / universe 13,552) ==="
FIG4A_KEEP_TIER34=1 python scripts/figures/fig4a_overview_candidates.py

echo "FIG4A_CAND_JOB_DONE $(date)"
