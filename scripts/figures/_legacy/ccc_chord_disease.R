#!/usr/bin/env Rscript
# ============================================================================
# ccc_chord_disease.R
#
# DISEASE-DYNAMIC companion to the static full-network CCC chord
# (ccc_full_network_chord.R -> fig3g_ccc_chord_all_lr.pdf).
#
# Same topology, recolored by disease direction:
#   * Ribbon WIDTH  = consensus L-R count per directed sender->receiver edge,
#                     computed IDENTICALLY to the static chord (specificity_rank
#                     <= SPEC_CUT present, consensus if present in >= DONOR_FRAC
#                     of donors, width = uniqueN(lr_pair) per (source,target)).
#   * Ribbon COLOR  = disease direction of that edge, aggregated across its L-R
#                     pairs from the stage LMM disease-beta (term ==
#                     "disease_stage_coarseSteatohepatitis", Estimate column):
#                        RED  = disease-gained (net positive beta)
#                        BLUE = disease-lost   (net negative beta)
#                        GRAY = no significant net change (q_within_ct >= 0.05
#                               for all LR pairs on the edge).
#
# Width encodes the communication topology; color encodes the disease change.
# The user places this beside fig3g_ccc_chord_all_lr.pdf in Illustrator; this
# script writes a STANDALONE panel only (no combined layout).
#
# Output (figures/main/fig3_RNAseq/panels/):
#   fig3g_ccc_chord_disease.pdf
#
# Run: micromamba activate rnaseq && Rscript scripts/figures/ccc_chord_disease.R
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
V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
META_V2 <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv")

PANEL_DIR <- file.path(FIG2_DIR, "panels")   # FIG2_DIR == .../fig3_RNAseq
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_DIS <- file.path(PANEL_DIR, "fig3g_ccc_chord_disease.pdf")

# Tunable thresholds (must MATCH the static chord for identical topology)
SPEC_CUT   <- 0.05   # LIANA specificity_rank cutoff for "present in a donor"
DONOR_FRAC <- 0.25   # fraction of donors required for a "consensus" interaction
Q_SIG      <- 0.05   # q_within_ct significance cutoff for the disease-beta mask

# Disease-direction ribbon palette (diverging; all TEXT stays black elsewhere)
COL_GAIN <- "#C2185B"   # RED  — disease-gained (net positive disease-beta)
COL_LOST <- "#1565C0"   # BLUE — disease-lost   (net negative disease-beta)
COL_NS   <- "#9E9E9E"   # GRAY — no significant net change (control gray)

# short cell-type labels + colours, identical to the static chord
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
# (a) TOPOLOGY  — consensus edge widths (identical to the static chord)
# ============================================================================
cat("\n[disease-chord] reading full donor LR matrix (topology)\n")
lr <- fread(file.path(V2_DIR, "all_donor_lr_scores_v2.tsv.gz"))

# protocol-contamination remediation: drop exclude_stage_analysis donors
meta <- fread(META_V2)
if (!"exclude_stage_analysis" %in% names(meta)) meta[, exclude_stage_analysis := FALSE]
keep_samples <- meta[exclude_stage_analysis != TRUE, unique(sample)]
lr <- lr[sample %in% keep_samples]
n_donors <- uniqueN(lr$sample)
cat(sprintf("  donors retained: %d\n", n_donors))

lr[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr[, present := specificity_rank <= SPEC_CUT]

cons <- lr[present == TRUE,
           .(n_donors_present = uniqueN(sample)),
           by = .(source, target, lr_pair)]
cons[, frac := n_donors_present / n_donors]
cons <- cons[frac >= DONOR_FRAC]
cat(sprintf("  consensus LR interactions (>=%.0f%% donors, spec<=%.2f): %d\n",
            100 * DONOR_FRAC, SPEC_CUT, nrow(cons)))

# ribbon width = number of consensus LR pairs per directed ct-pair
adj <- cons[, .(value = uniqueN(lr_pair)), by = .(source, target)]

# ============================================================================
# (b) DISEASE DIRECTION  — aggregate the stage-LMM disease-beta per edge
# ============================================================================
cat("[disease-chord] reading stage LMM disease-beta (v3 coarse)\n")
lmm <- fread(file.path(V3_DIR, "stage_lr_lmm_coarse_v3.tsv"))
lmm <- lmm[term == "disease_stage_coarseSteatohepatitis"]
cat(sprintf("  disease-beta rows (Steatohepatitis term): %d\n", nrow(lmm)))

# significance mask: an LR pair counts as a "significant" disease change only
# when q_within_ct < Q_SIG; non-significant pairs contribute 0 to the net.
lmm[, sig := !is.na(q_within_ct) & q_within_ct < Q_SIG]
lmm[, dir := fifelse(sig & Estimate > 0,  1L,
              fifelse(sig & Estimate < 0, -1L, 0L))]

# per directed edge: mean disease-beta (magnitude) + net signed count of
# significant gained/lost LR pairs (used for the RED / BLUE / GRAY call).
edge_dir <- lmm[, .(
  mean_beta      = mean(Estimate, na.rm = TRUE),
  n_lr           = .N,
  n_sig          = sum(sig),
  n_sig_gained   = sum(dir ==  1L),
  n_sig_lost     = sum(dir == -1L),
  net_signed     = sum(dir),                       # (#gained - #lost)
  mean_sig_beta  = { b <- Estimate[sig]; if (length(b)) mean(b) else NA_real_ }
), by = .(source, target)]

# join disease direction onto the topology edges
adj <- merge(adj, edge_dir, by = c("source", "target"), all.x = TRUE)

# direction call: net signed count of significant LR pairs decides color;
# ties / all-non-significant edges -> GRAY.
adj[, direction := fifelse(!is.na(net_signed) & net_signed > 0, "gained",
                    fifelse(!is.na(net_signed) & net_signed < 0, "lost",
                            "ns"))]
adj[, ribbon_col := fifelse(direction == "gained", COL_GAIN,
                     fifelse(direction == "lost",  COL_LOST, COL_NS))]

# short labels for circlize
adj[, from := ct_short_map[source]]
adj[, to   := ct_short_map[target]]
adj <- adj[!is.na(from) & !is.na(to)]
setorder(adj, -value)

n_g <- adj[direction == "gained", .N]
n_l <- adj[direction == "lost",   .N]
n_n <- adj[direction == "ns",     .N]
cat(sprintf("  edges: %d total | %d disease-gained (RED) | %d disease-lost (BLUE) | %d non-sig (GRAY)\n",
            nrow(adj), n_g, n_l, n_n))

# persist the edge table for provenance
fwrite(adj[, .(source, target, from, to, n_consensus_lr = value,
               mean_beta, n_lr, n_sig, n_sig_gained, n_sig_lost,
               net_signed, mean_sig_beta, direction, ribbon_col)],
       file.path(DATA_DIR, "ccc_chord_disease_data.csv"))

# report the strongest disease gain / loss edges (by net signed count, then beta)
strongest_gain <- adj[direction == "gained"][order(-net_signed, -mean_sig_beta)]
strongest_loss <- adj[direction == "lost"  ][order(net_signed,  mean_sig_beta)]
cat("\n  Strongest DISEASE-GAINED edges (sender -> receiver):\n")
if (nrow(strongest_gain)) {
  print(head(strongest_gain[, .(edge = paste(from, "->", to),
                                n_consensus_lr = value, net_signed,
                                mean_sig_beta = round(mean_sig_beta, 3))], 5))
} else cat("    (none)\n")
cat("\n  Strongest DISEASE-LOST edges (sender -> receiver):\n")
if (nrow(strongest_loss)) {
  print(head(strongest_loss[, .(edge = paste(from, "->", to),
                                n_consensus_lr = value, net_signed,
                                mean_sig_beta = round(mean_sig_beta, 3))], 5))
} else cat("    (none)\n")

# ============================================================================
# (c) RENDER  — reuse pretty_chord() with the static chord's node order/style
# ============================================================================
# sector order: senders by descending out-degree (identical to static chord)
node_order <- ct_short_map[c("Hepatocytes","Endothelial cells","Fibroblasts",
                             "Macrophages","Cholangiocytes","T cells",
                             "B cells","Resident NK","Plasma cells")]
nodes <- node_order[node_order %in% c(adj$from, adj$to)]
grid_col <- setNames(ct_short_col[nodes], nodes)

# ribbon colour = disease direction (semi-transparent to match static alpha 0.62)
adj[, col := alpha(ribbon_col, 0.62)]

FIG_W <- 95 / 25.4
FIG_H <- 96 / 25.4
cairo_pdf(OUT_DIS, width = FIG_W, height = FIG_H)
layout(matrix(c(1, 2), nrow = 2), heights = c(86, 10) / 25.4)
par(mar = c(0.5, 0.5, 0.5, 0.5), family = "Helvetica")

pretty_chord(adj, nodes, grid_col,
             start.degree = 90, gap.degree = 4,
             show_axis = FALSE, show_total = TRUE,
             label_fontsize = 6, total_fontsize = 6,
             text_family = "Helvetica", label_font = 1, total_font = 1,
             label_offset = 5.5, total_offset = 1.4, canvas_lim = 1.28)

par(mar = c(0, 0, 0, 0), family = "Helvetica"); plot.new()
legend("center",
       legend = c("Disease-gained (β > 0)",
                  "Disease-lost (β < 0)",
                  "No significant change"),
       fill = c(COL_GAIN, COL_LOST, COL_NS),
       border = NA, horiz = TRUE, bty = "n", cex = 0.5, x.intersp = 0.4,
       text.col = "black",
       title = "Ribbon colour = disease direction", title.font = 1)
dev.off()
cat(sprintf("\n[disease-chord] wrote %s\n", OUT_DIS))

# caption (house style: descriptive text via message(), never a subtitle)
message(
  "CAPTION fig3g_ccc_chord_disease: Disease-dynamic cell-cell communication ",
  "chord. Ribbon WIDTH = number of consensus ligand-receptor interactions per ",
  "directed sender->receiver cell-type edge (LIANA specificity_rank <= 0.05, ",
  "consensus in >= 25% of donors) - the identical topology to the static ",
  "full-network chord. Ribbon COLOUR = disease direction from the stage mixed ",
  "model (term = disease_stage_coarseSteatohepatitis): RED = disease-gained ",
  "(net positive disease-beta among q_within_ct < 0.05 LR pairs), BLUE = ",
  "disease-lost (net negative), GREY = no significant net change. Sector ",
  "colours and node order match the static chord.")

cat("\nDONE\n")
