#!/bin/bash
#SBATCH --job-name=liver_rerun_02_12
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/rerun_02_12_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/rerun_02_12_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=6:00:00

# Full re-run of integration pipeline (scripts 02-12).
# Generic cascade wrapper for after metadata/yaml or normalization changes
# (e.g., the 2026-05-15 PRJNA512027 removal).

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== Pipeline re-run started: $(date) ==="
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK | Mem: 400G | Partition: bigmem"
echo ""

run_step() {
  local num="$1"
  local name="$2"
  echo ""
  echo "============================================================"
  echo "  Step $num: $name"
  echo "============================================================"
  Rscript "${SCRIPTS}/${name}" 2>&1
  if [ $? -ne 0 ]; then
    echo "FAILED: $name at $(date)"
    exit 1
  fi
  echo "  Completed: $(date)"
}

run_step "02" "02_per_study_de.R"
run_step "03" "03_integrate_counts.R"
run_step "04" "04_variance_partition.R"
run_step "05" "05_dream_mega_analysis.R"
run_step "07" "07_consensus_degs.R"
run_step "08" "08_pathway_analysis.R"
run_step "09" "09_batch_correction_umap.R"
run_step "10" "10_volcano_plots.R"
run_step "12" "12_library_intersection.R"

echo ""
echo "============================================================"
echo "  ALL STEPS COMPLETE: $(date)"
echo "============================================================"
