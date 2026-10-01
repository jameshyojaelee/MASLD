#!/bin/bash
#SBATCH --job-name=fig4e_corr
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4e_corr_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4e_corr_%j.err

# Figure 4E on the corrected-control stage release (review item 1, Phase 3).
#
# The transport table is not written by build_pi_stage_extensions.R. Its chain is
#   render_pi_singlecell_panels.R (STAGE_RELEASE_ROOT)   -> fig4e_bulk_projection.tsv
#   05_build_figure_sources.py (WORKED_LOCI_BULK_PROJECTION) -> fig4f_bulk_stage_transport.tsv
#   render_fig4def_multicellular_story.R:146-160         -> annotation columns
#   render_fig4e_bulk_uniform_points.R                   -> panel
# run_worked_loci_corrected.sbatch already runs the first two steps on the corrected
# stage directory. This job annotates its table as render_fig4def_multicellular_story.R
# does (with the BH family size read from the stage release, not the literal
# 23,370), renders 4E, and compares all 8 points with the adopted panel.
#
#   sbatch --dependency=afterok:<worked-loci job> --export=ALL,\
# WORKED_LOCI_ROOT=<the worked-loci CANDIDATE_ROOT>,\
# BULK_STAGE_DIR=<bulk release>/analysis/stage_extensions,\
# FIG4E_ROOT=<BASE>/figures/candidates/fig4e-corrected-controls-<UTC ts> \
#     scripts/figures/run_fig4e_corrected_controls.sbatch
set -euo pipefail
umask 027
: "${WORKED_LOCI_ROOT:?Set the worked-loci candidate root built on the corrected stage release}"
: "${BULK_STAGE_DIR:?Set the corrected stage_extensions directory used by that worked-loci run}"
: "${FIG4E_ROOT:?Set a new figures/candidates/fig4e-corrected-controls-<UTC ts> directory}"
test ! -e "$FIG4E_ROOT"
project_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
rnaseq_env=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export PATH="/gpfs/commons/home/jameslee/.local/bin:${PATH}"
export PYTHONNOUSERSITE=1
transport="$WORKED_LOCI_ROOT/source_tables/fig4f_bulk_stage_transport.tsv"
adopted="$project_root/figures/main/fig4_singlecell_programs/source_tables/current_candidate/fig4e_bulk_tissue_state_transport.tsv"
test -s "$transport"
test -s "$WORKED_LOCI_ROOT/upstream/fig4e_projection_corrected/source_tables/fig4e_bulk_projection.tsv"
test -s "$BULK_STAGE_DIR/stage_extension_all_gene_results.tsv"
mkdir -p "$FIG4E_ROOT/inputs"

"$rnaseq_env/bin/Rscript" --vanilla - "$transport" "$BULK_STAGE_DIR/stage_extension_all_gene_results.tsv" \
  "$FIG4E_ROOT/inputs/fig4e_bulk_tissue_state_transport.tsv" <<'RS'
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
f <- fread(args[[1]], na.strings = c("", "NA"))
stage <- fread(args[[2]], select = c("contrast", "bh_family_size"))
family <- unique(stage[contrast %in% unique(f$contrast), bh_family_size])
stopifnot(nrow(f) == 8L, setequal(f$contrast, paste0("F", 1:4, "_vs_F0")), length(family) == 1L)
# Same annotation as render_fig4def_multicellular_story.R:146-160.
f[, biological_unit := "bulk RNA-seq participant"]
f[, multiple_testing_family := paste0(
  format(family, big.mark = ","),
  " genes within each fibrosis-stage-versus-F0 contrast; no program-level significance test"
)]
f[, evidence_state := "descriptive_tissue_transport"]
f[, state_reason := paste0(
  "fixed L1-weighted member-gene effect with ", n_bh, "/", n_genes,
  " available members passing the gene-level BH family"
)]
f[, unresolved_alternative := paste0(
  "bulk tissue composition and shared fibrosis state cannot establish cellular origin"
)]
f[, claim_boundary := paste0(
  "cross-sectional fixed-member tissue-state transport; not longitudinal progression or a program-level inferential test"
)]
fwrite(f, args[[3]], sep = "\t")
RS

FIG4E_TRANSPORT_SOURCE="$FIG4E_ROOT/inputs/fig4e_bulk_tissue_state_transport.tsv" \
FIG4E_UNIFORM_OUT="$FIG4E_ROOT/render" \
  "$rnaseq_env/bin/Rscript" --vanilla "$project_root/scripts/figures/render_fig4e_bulk_uniform_points.R"
pdf="$FIG4E_ROOT/render/panels/fig4e_bulk_tissue_state_transport.pdf"
test "$(pdfinfo "$pdf" | awk '/^Pages:/ {print $2}')" = "1"
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "$pdf"
mkdir -p "$FIG4E_ROOT/proofs"
pdftotext "$pdf" "$FIG4E_ROOT/proofs/fig4e_bulk_tissue_state_transport.txt"
pdftoppm -png -r 400 -singlefile "$pdf" "$FIG4E_ROOT/proofs/fig4e_bulk_tissue_state_transport"

# All 8 plotted points, adopted versus corrected.
"$rnaseq_env/bin/python3" - "$adopted" "$FIG4E_ROOT/render/source_tables/fig4e_bulk_tissue_state_transport.tsv" \
  "$FIG4E_ROOT/fig4e_old_vs_new.tsv" <<'PY'
import sys
import pandas as pd
keys = ["program_uid", "stage"]
cols = ["effect", "n_genes", "n_bh", "bh_weight_fraction", "direction_agreement"]
old, new = (pd.read_csv(p, sep="\t")[keys + ["program_name"] + cols] for p in sys.argv[1:3])
m = old.merge(new.drop(columns="program_name"), on=keys, suffixes=("_adopted", "_corrected"),
              how="outer", validate="one_to_one")
assert len(m) == 8 and m.notna().all().all(), "the 8 program-by-stage points do not match"
m["effect_delta"] = m["effect_corrected"] - m["effect_adopted"]
m.sort_values(keys).to_csv(sys.argv[3], sep="\t", index=False)
print(m[["program_name", "stage", "effect_adopted", "effect_corrected", "effect_delta",
         "n_genes_adopted", "n_genes_corrected", "n_bh_adopted", "n_bh_corrected"]].to_string(index=False))
PY

find "$FIG4E_ROOT" -type f ! -name input_output_manifest.tsv -print0 | sort -z | \
  xargs -0 sha256sum > "$FIG4E_ROOT/input_output_manifest.tsv"
sha256sum "$transport" "$adopted" >> "$FIG4E_ROOT/input_output_manifest.tsv"
echo "FIG4E_CORRECTED_COMPLETE $FIG4E_ROOT"
