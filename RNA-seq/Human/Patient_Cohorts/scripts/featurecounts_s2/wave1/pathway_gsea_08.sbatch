#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/pathway_gsea_08_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/pathway_gsea_08_%j.err
# Wave1 [H] ATLAS FEEDER (CORRECTED): 08_pathway_analysis.R writes results/integration/gsea_results.csv,
# which is what 27a:1057 reads. (User named 08b, but 08b writes gsea_kallisto/gsea_all_contrasts.csv,
# NOT consumed by 27a -- see notes.) Reads canonical -s2 dream_results.csv. Env=rnaseq (fgsea/clusterProfiler).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration
echo "=== 08_pathway_analysis (gsea_results.csv, atlas feeder) | $(date) ==="
Rscript scripts/08_pathway_analysis.R
if [ $? -ne 0 ]; then echo "FAILED: 08_pathway_analysis"; exit 1; fi
echo "=== done $(date) ==="
