#!/bin/bash
#SBATCH --job-name=workedloci
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/aging_inspired_resource_figures/slurm/worked_loci_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/aging_inspired_resource_figures/slurm/worked_loci_%j.log

# Figure 4-6 worked-loci rebuild against the corrected bulk and COLOC releases.
# Required: CANDIDATE_ROOT (must not exist), COLOC_RELEASE_DIR, BULK_STAGE_DIR,
# PROTEIN_DATA_DIR. The GNMT genetics node is computed from the COLOC tier-1/2
# table; it is not hand-set.
#   COLOC_RELEASE_DIR  susie_coloc_r2_<ts>/ (gene_level_coloc.csv, susie_coloc_all_gwas.csv,
#                      tier12/gene_level_coloc_tier12.csv)
#   BULK_STAGE_DIR     <bulk release>/analysis/stage_extensions (stage_all_gene_results.tsv,
#                      cohort_disease_all_gene_results.tsv, stage_extension_all_gene_results.tsv)
#   PROTEIN_DATA_DIR   <fig5 candidate>/figure5/data (composite_mrna_protein_corrected.csv,
#                      composite_protein_histology_partial_corrected.csv)
# Optional: COLOC_TIER12 and COLOC_MASTER override the r2 paths.
#
# The Figure 4F bulk transport (fig4e_bulk_projection.tsv) is recomputed here with
# render_pi_singlecell_panels.R on BULK_STAGE_DIR. Its F-versus-F0 fits select by
# fibrosis stage, not by control label, but the corrected refit re-filters the gene
# universe with filterByExpr(~ 0 + group_binary). The producer is first rerun on
# the stage release behind the retained fig4 v4/v5 table and must reproduce it.

set -euo pipefail
umask 027

: "${CANDIDATE_ROOT:?CANDIDATE_ROOT is required}"
: "${COLOC_RELEASE_DIR:?COLOC_RELEASE_DIR is required}"
: "${BULK_STAGE_DIR:?BULK_STAGE_DIR is required}"
: "${PROTEIN_DATA_DIR:?PROTEIN_DATA_DIR is required}"
if [[ -e "$CANDIDATE_ROOT" ]]; then
  echo "Refusing existing CANDIDATE_ROOT: $CANDIDATE_ROOT" >&2
  exit 2
fi

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export FIG_CAND_ROOT="$CANDIDATE_ROOT"
export WORKED_LOCI_COLOC_GENE_LEVEL="$COLOC_RELEASE_DIR/gene_level_coloc.csv"
export WORKED_LOCI_COLOC_TIER12="${COLOC_TIER12:-$COLOC_RELEASE_DIR/tier12/gene_level_coloc_tier12.csv}"
export WORKED_LOCI_COLOC_MASTER="${COLOC_MASTER:-$COLOC_RELEASE_DIR/susie_coloc_all_gwas.csv}"
export WORKED_LOCI_STAGE_DIR="$BULK_STAGE_DIR"
export WORKED_LOCI_PROTEIN_DATA_DIR="$PROTEIN_DATA_DIR"
for f in "$WORKED_LOCI_COLOC_GENE_LEVEL" "$WORKED_LOCI_COLOC_TIER12" "$WORKED_LOCI_COLOC_MASTER" \
         "$BULK_STAGE_DIR/stage_all_gene_results.tsv" \
         "$BULK_STAGE_DIR/stage_extension_all_gene_results.tsv" \
         "$BULK_STAGE_DIR/cohort_disease_all_gene_results.tsv" \
         "$PROTEIN_DATA_DIR/composite_mrna_protein_corrected.csv" \
         "$PROTEIN_DATA_DIR/composite_protein_histology_partial_corrected.csv"; do
  test -s "$f" || { echo "Missing input: $f" >&2; exit 2; }
done

export SCRIPT_ROOT=${MASLD_PROJECT_ROOT}/scripts/manuscript/aging_inspired_resource_figures
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export PYTHONNOUSERSITE=1
python_bin=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python

mkdir "$CANDIDATE_ROOT"
echo "[worked-loci] host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -u)"
echo "[worked-loci] tier12=$WORKED_LOCI_COLOC_TIER12 stage=$BULK_STAGE_DIR protein=$PROTEIN_DATA_DIR"

# Figure 4F bulk projection: reproduce the retained table, then recompute it on the
# corrected stage release. Inputs match the fig4 v4 render (its provenance sbatch).
rscript=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
fig4_v5=${MASLD_PROJECT_ROOT}/figures/candidates/fig4-stage-terminology-corrected-2026-08-13-v5
project_fig4e() {  # <stage dir> <output root>
  mkdir -p "$2"
  FIGURE_CANDIDATE_ROOT="$2" STAGE_RELEASE_ROOT="$1" \
  PROGRAM_RELEASE_ROOT=${MASLD_PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot \
  CCC_TRAJECTORY_FILE="$fig4_v5/source_tables/fig4d_communication_trajectories.tsv" \
  CCC_AUDIT_FILE="$fig4_v5/source_tables/fig4d_communication_audit.tsv" \
    "$rscript" --vanilla "${MASLD_PROJECT_ROOT}/scripts/figures/render_pi_singlecell_panels.R"
}
project_fig4e "${MASLD_PROJECT_ROOT}/figures/candidates/pi-figure-redesign-2026-08-13-v6/analysis/stage_extensions" \
  "$CANDIDATE_ROOT/upstream/fig4e_projection_reproduction"
"${python_bin}" - "$fig4_v5/source_tables/fig4e_bulk_projection.tsv" \
  "$CANDIDATE_ROOT/upstream/fig4e_projection_reproduction/source_tables/fig4e_bulk_projection.tsv" <<'PY'
import sys
import numpy as np
import pandas as pd
keys = ["program_uid", "contrast"]
old, new = (pd.read_csv(p, sep="\t").sort_values(keys).reset_index(drop=True) for p in sys.argv[1:3])
assert old[keys].equals(new[keys]), "projection rows differ from the retained table"
for c in ("effect", "n_genes", "n_bh", "bh_weight_fraction", "direction_agreement"):
    assert np.allclose(old[c], new[c], rtol=1e-10, atol=1e-12), f"projection does not reproduce: {c}"
print(f"[worked-loci] fig4e projection reproduced: {len(new)} rows")
PY
project_fig4e "$BULK_STAGE_DIR" "$CANDIDATE_ROOT/upstream/fig4e_projection_corrected"
export WORKED_LOCI_BULK_PROJECTION="$CANDIDATE_ROOT/upstream/fig4e_projection_corrected/source_tables/fig4e_bulk_projection.tsv"
"${python_bin}" "${SCRIPT_ROOT}/05_build_figure_sources.py"
# 08 is the compact renderer that 07 validates; it also writes sessionInfo.txt.
micromamba run -n decontx Rscript "${SCRIPT_ROOT}/08_render_compact_panels.R"
"${python_bin}" -m pip freeze > "${CANDIDATE_ROOT}/pip_freeze.txt"
"${python_bin}" "${SCRIPT_ROOT}/07_validate_figure_candidate.py"
echo "[worked-loci] done=$(date -u)"
