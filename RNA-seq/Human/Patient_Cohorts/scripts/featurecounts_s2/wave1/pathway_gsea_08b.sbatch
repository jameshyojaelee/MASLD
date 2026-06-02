#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/pathway_gsea_08b_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/pathway_gsea_08b_%j.err
# Wave1 [H] COMPREHENSIVE (user-named): 08b_unified_pathway_gsea.R -> gsea_kallisto/gsea_all_contrasts.csv
# (multi-contrast x multi-collection). NOT the 27a feeder (that is 08 -> gsea_results.csv). Reads -s2
# dream_results.csv. Env=rnaseq (per its own run_08b_unified_gsea.sh; NOT the isolated pathway env).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration
echo "=== 08b_unified_pathway_gsea (comprehensive) | $(date) ==="
Rscript scripts/08b_unified_pathway_gsea.R
if [ $? -ne 0 ]; then echo "FAILED: 08b_unified_pathway_gsea"; exit 1; fi
echo "=== done $(date) ==="
