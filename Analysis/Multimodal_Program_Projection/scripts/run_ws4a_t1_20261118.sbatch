#!/usr/bin/env bash
#SBATCH --job-name=ws4a_t1_prot
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/ws4a_t1_prot_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/ws4a_t1_prot_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

export PROTEOMICS_OUT_DIR="$PWD/Analysis/Multimodal_Program_Projection/results/proteomics_20261118"
mkdir -p "$PROTEOMICS_OUT_DIR"

# 1. Core proteomics projection (58-sample DE + 46-MASLD-only v1 22-program & histology)
"$R_BIN" --vanilla Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection_20261118.R

# 2. Unified v1 (22) + v2 (117) MASLD histology burden with Dirichlet(1,...,1) weights
"$R_BIN" --vanilla figures/candidates/panel-deep-fixes-20261017T035633Z/fig5/code/fig5e_protein_histology_117_20261118.R

# 3. Trio covariance (GNMT-MAT1A-CYP2C19) with 46-MASLD qnorm + Dirichlet bootstrap
"$R_BIN" --vanilla scripts/analysis/fig34_corrections/02_protein_20261118.R "$PROTEOMICS_OUT_DIR/trio_covariance"

# 4. Fig 5B protein triage verification written cleanly under proteomics_20261118/figure5/panels/data
export FIG5B_PROTEIN_DE="$PROTEOMICS_OUT_DIR/protein_de_adjusted.tsv"
export FIG5B_BULK_DEG="$PWD/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
export FIGURE_CANDIDATE_ROOT="$PROTEOMICS_OUT_DIR"
"$R_BIN" --vanilla scripts/figures/fig5b_protein_triage.R
