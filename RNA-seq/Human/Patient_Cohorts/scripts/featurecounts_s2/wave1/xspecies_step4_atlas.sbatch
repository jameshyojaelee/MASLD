#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_step4_atlas_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_step4_atlas_%j.err
# Wave1 [E] step4 (ATLAS FEEDER): 04 viz -> 05 unified atlas -> concordance_atlas_unified.csv (27a).
# Runs after step2 (02 array) AND step3 (03 network) complete.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
S=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== xspecies 04 visualizations | $(date) ==="
Rscript "$S/04_concordance_visualizations.R"
if [ $? -ne 0 ]; then echo "FAILED: 04_concordance_visualizations"; exit 1; fi
echo "=== xspecies 05 unified concordance atlas (FEEDER) | $(date) ==="
Rscript "$S/05_unified_concordance_atlas.R"
if [ $? -ne 0 ]; then echo "FAILED: 05_unified_concordance_atlas"; exit 1; fi
echo "=== done $(date) ==="
