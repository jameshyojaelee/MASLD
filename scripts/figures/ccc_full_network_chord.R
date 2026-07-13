#!/usr/bin/env Rscript
# ============================================================================
# ccc_full_network_chord.R
#
# Full LIANA communication-network chord for the single-cell CCC analysis.
# Does NOT touch the canonical ccc_v3_panels.R outputs (the stage-gated chord).
#
# Output (figures/main/fig3_RNAseq/panels/):
#   figs3_ccc_chord_all_lr.pdf   (supp; was fig3g, demoted 2026-07-02)
#      FULL LIANA communication network. Every directed sender->receiver
#      cell-type pair across all 9 cell types. Ribbon width = number of
#      *consensus* ligand-receptor interactions, where an LR pair is "present"
#      in a donor when specificity_rank <= SPEC_CUT and "consensus" when present
#      in >= DONOR_FRAC of the donors. Ribbon colored by SENDER cell type.
#      This is the unfiltered communication landscape (NOT the stage-gated set).
#
# Run: ~/micromamba/envs/rnaseq/bin/Rscript scripts/figures/ccc_full_network_chord.R
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/ccc_chord_style.R"))

V2_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2")
META_V2 <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv")

PANEL_DIR <- file.path(FIG2_DIR, "panels")   # FIG2_DIR == .../fig3_RNAseq
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_ALL  <- file.path(PANEL_DIR, "figs3_ccc_chord_all_lr.pdf")   # DEMOTED to supp 2026-07-02 (was fig3g)

# Tunable thresholds (printed to log for reproducibility)
SPEC_CUT   <- 0.05   # LIANA specificity_rank cutoff for "present in a donor"
DONOR_FRAC <- 0.25   # fraction of donors required for a "consensus" interaction

ct_short_map <- c(
  "Hepatocytes"       = "Hep",
  "Cholangiocytes"    = "Chol",
  "Endothelial cells" = "Endo",
  "Fibroblasts"       = "Fib",
  "Macrophages"       = "Mac",
  "T cells"           = "Tcell",
  "B cells"           = "Bcell",
  "Resident NK"       = "NK",
  "Plasma cells"      = "Plasma"
)
ct_short_col <- setNames(ct_palette[names(ct_short_map)], ct_short_map)

# ============================================================================
# PANEL 1: ccc_chord_all_lr.pdf  — full consensus communication network
# ============================================================================
cat("\n[all-LR] reading full donor LR matrix\n")
lr <- fread(file.path(V2_DIR, "all_donor_lr_scores_v2.tsv.gz"))

# protocol-contamination remediation: drop exclude_stage_analysis donors
meta <- fread(META_V2)
if (!"exclude_stage_analysis" %in% names(meta)) meta[, exclude_stage_analysis := FALSE]
keep_samples <- meta[exclude_stage_analysis != TRUE, unique(sample)]
lr <- lr[sample %in% keep_samples]
n_donors <- uniqueN(lr$sample)   # donors actually present in the LR matrix
cat(sprintf("  donors retained: %d\n", n_donors))
cat(sprintf("  raw donor x LR x ct-pair rows: %d\n", nrow(lr)))

lr[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr[, present := specificity_rank <= SPEC_CUT]

# consensus interaction = LR pair present in >= DONOR_FRAC of donors for a ct-pair
cons <- lr[present == TRUE,
           .(n_donors_present = uniqueN(sample)),
           by = .(source, target, lr_pair)]
cons[, frac := n_donors_present / n_donors]
cons <- cons[frac >= DONOR_FRAC]
cat(sprintf("  consensus LR interactions (>=%.0f%% donors, spec<=%.2f): %d\n",
            100 * DONOR_FRAC, SPEC_CUT, nrow(cons)))

# ribbon width = number of consensus LR pairs per directed ct-pair
adj_all <- cons[, .(value = uniqueN(lr_pair)), by = .(source, target)]
adj_all[, from := ct_short_map[source]]
adj_all[, to   := ct_short_map[target]]
adj_all <- adj_all[!is.na(from) & !is.na(to)]
setorder(adj_all, -value)

fwrite(adj_all[, .(source, target, from, to, n_consensus_lr = value)],
       file.path(DATA_DIR, "ccc_chord_all_lr_data.csv"))
fwrite(cons[order(source, target, -frac)],
       file.path(DATA_DIR, "ccc_chord_all_lr_consensus_pairs.csv"))

# sector order: senders by descending out-degree for a clean visual flow
node_order <- ct_short_map[c("Hepatocytes","Endothelial cells","Fibroblasts",
                             "Macrophages","Cholangiocytes","T cells",
                             "B cells","Resident NK","Plasma cells")]
nodes_all <- node_order[node_order %in% c(adj_all$from, adj_all$to)]
grid_col_all <- setNames(ct_short_col[nodes_all], nodes_all)

# ribbon color = sender cell type (semi-transparent)
adj_all[, col := alpha(ct_short_col[from], 0.62)]

FIG_W <- 95 / 25.4
FIG_H <- 96 / 25.4
cairo_pdf(OUT_ALL, width = FIG_W, height = FIG_H)
layout(matrix(c(1, 2), nrow = 2), heights = c(86, 10) / 25.4)
par(mar = c(0.5, 0.5, 0.5, 0.5), family = "Helvetica")

pretty_chord(adj_all, nodes_all, grid_col_all,
             start.degree = 90, gap.degree = 4,
             show_axis = FALSE, show_total = TRUE,
             label_fontsize = 6, total_fontsize = 6,
             text_family = "Helvetica", label_font = 1, total_font = 1,
             label_offset = 5.5, total_offset = 1.4, canvas_lim = 1.28)

par(mar = c(0, 0, 0, 0), family = "Helvetica"); plot.new()
legend("center", legend = names(grid_col_all), fill = grid_col_all,
       border = NA, horiz = TRUE, bty = "n", cex = 0.5, x.intersp = 0.4,
       title = "Sender / receiver cell type", title.font = 1)
dev.off()
cat(sprintf("[all-LR] wrote %s\n", OUT_ALL))

cat("\nDONE\n")
