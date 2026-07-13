#!/bin/bash
#SBATCH --job-name=heatmap
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/extreme_hm_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/extreme_hm_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=0-6

# Extreme-phenotype heatmap N-sweep -- one DEG count per array task, fanned
# across nodes. n=200 writes the default (no-suffix) file; the rest write
# _nN variants for pick-and-choose. cpu partition, default QOS.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u; micromamba activate rnaseq; set -u

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$ROOT"

NS=(50 100 150 200 300 500 1000)
N=${NS[$SLURM_ARRAY_TASK_ID]}
if [ "$N" = "200" ]; then SUF=""; else SUF="_n${N}"; fi
echo "=== extreme heatmap n=$N (suffix '$SUF') | job $SLURM_ARRAY_JOB_ID task $SLURM_ARRAY_TASK_ID | $(date) ==="

TOP_N_DEG=$N HEATMAP_OUT_SUFFIX="$SUF" \
  Rscript scripts/figures/figS_extreme_phenotype_heatmap.R 2>&1

echo "=== done n=$N: $(date) ==="
