#!/usr/bin/env bash
#SBATCH --job-name=A10b_strat_gsea
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=outputs/team_A/A10b_launch_slurm_%j.out
#SBATCH --error=outputs/team_A/A10b_launch_slurm_%j.err

# activate environment BEFORE set -eu — conda activate scripts reference
# unset variables (ADDR2LINE etc.) which trip nounset.
# shellcheck source=/dev/null
source ~/.bashrc
micromamba activate rnaseq

set -eu

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "[$(date +%H:%M:%S)] A10b stratified GSEA launch"
echo "Host: $(hostname)"
echo "PWD : $(pwd)"
echo "Job : ${SLURM_JOB_ID:-local}"

Rscript RNA-seq/218c2_pathway_gsea_stratified.R

echo "[$(date +%H:%M:%S)] A10b stratified GSEA done"
