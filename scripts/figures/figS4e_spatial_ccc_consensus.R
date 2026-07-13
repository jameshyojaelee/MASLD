#!/usr/bin/env Rscript
# ==============================================================================
# Fig S4e (DEMOTED from main 4h 2026-07-08) — spatial grounding of inferred cell-cell signaling (Visium spot-adjacency,
# GSE192741). Ligand-receptor pairs inferred from DISSOCIATED scRNA-seq (LIANA;
# Fig 3's RNA-only lens) that are ALSO detected in INTACT TISSUE by a
# spatial-proximity method (squidpy on Visium spot adjacency). SET-LEVEL over the
# pairs found by BOTH; disease-DIRECTION concordance is defined because GSE192741
# carries a healthy comparator.
#   p_sc  = LIANA disease-score (x) vs squidpy Δ (y), one point per pair, coloured
#           by direction-concordance (same disease-sign in both methods).
#   p_bar = source-attribution skew: which cell type each method calls the SOURCE
#           (squidpy → hepatocyte-autocrine-leaning; LIANA → non-parenchymal/immune).
#
# PAIR-LEVEL recovery of inferred signaling, NOT cell-type-resolved proof. Vu et al.
# Visium is deliberately NOT included here (no healthy comparator, no per-spot
# fibrosis staging → cannot support a disease-DIRECTION claim; the established
# non-pooling rule bars using GSE192741-healthy as its reference).
#
# The companion CosMx single-cell-proximity exhibit (Govaere 2026, IL32-hepatocyte
# -> macrophage) was split out to figS4c_il32_macrophage_cosmx.R 2026-07-08 (user
# call: two genuinely distinct exhibits crammed into one callout; the CosMx one is a
# narrower single-mechanism deep-dive better suited to supplement than main text).
#
# Data: Analysis/Spatial/results/spatial_ccc_multimethod/scRNA_spatial_LR_consensus.csv
# Output: figures/main/fig4_validation/panels/figS4e.pdf (supp; kept in panels/ per figS4 convention)
# Env: rnaseq
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ==============================================================================
# PART 1 — GSE192741 Visium: LIANA (scRNA) vs squidpy (spatial) direction-concordance
# ==============================================================================
d <- fread(file.path(BASE, "Analysis/Spatial/results/spatial_ccc_multimethod/scRNA_spatial_LR_consensus.csv"))
n_pairs <- nrow(d)
n_conc  <- sum(d$both_concordant)
pct     <- round(100 * n_conc / n_pairs, 1)
d[, conc := factor(ifelse(both_concordant, "direction-concordant", "discordant"),
                   levels = c("direction-concordant", "discordant"))]
d[, pair := paste0(ligand, "–", receptor)]

# source-attribution classes (which cell type each method calls the SOURCE)
hep_class <- function(x) ifelse(x == "Hepatocytes", "Hepatocyte", "Non-parenchymal")
liana_hep <- sum(d$liana_source == "Hepatocytes")
sq_hep    <- sum(d$sq_source   == "Hepatocytes")
sq_auto   <- sum(d$sq_pair_type == "autocrine")
imm <- c("cDC1s","cDC2s","pDCs","Mig.cDCs","Macrophages","Mono+mono derived cells",
         "T cells","B cells","Plasma cells","Resident NK","Circulating NK/NKT","Neutrophils","Basophils")
liana_imm <- sum(d$liana_source %in% imm)
cat(sprintf("[cccS4/P1] %d LR pairs in BOTH methods; %d direction-concordant (%.1f%%)\n", n_pairs, n_conc, pct))
cat(sprintf("[cccS4/P1] squidpy source hep=%d/%d (autocrine=%d); LIANA source hep=%d/%d (immune=%d)\n",
            sq_hep, n_pairs, sq_auto, liana_hep, n_pairs, liana_imm))

conc_col <- c("direction-concordant" = "#00695C", "discordant" = "#9E9E9E")
dc <- d[both_concordant == TRUE]
dc[, comb := frank(abs(liana_score_diff)) + frank(abs(sq_delta_expr))]
lab <- dc[order(-comb)][1:6]
# LAMC1-DAG1 sits in the dense bottom-left point cluster; nudge its label below the dots for legibility.
lab[, nudge_y := ifelse(pair == "LAMC1–DAG1", -0.65, 0)]

p_sc <- ggplot(d, aes(liana_score_diff, sq_delta_expr)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_point(aes(fill = conc), shape = 21, size = 1.7, stroke = 0.3, colour = "white") +
  ggrepel::geom_text_repel(data = lab, aes(label = pair), size = GEOM_TEXT_6PT, fontface = "italic",
                           colour = "black", max.overlaps = 40, min.segment.length = 0,
                           nudge_y = lab$nudge_y,
                           segment.size = 0.2, segment.colour = "grey70", box.padding = 0.32, seed = 42) +
  scale_fill_manual(values = conc_col, name = NULL) +
  scale_x_continuous(limits = c(-1.05, 1.05), breaks = c(-1, 0, 1)) +
  labs(x = "LIANA scRNA disease score",
       y = "squidpy spatial Δ") +
  ggtitle(sprintf("%d/%d (%.0f%%) direction-concordant", n_conc, n_pairs, pct)) +
  theme_masld() +
  theme(plot.title = element_text(size = 6, hjust = 0, face = "plain", lineheight = 0.95),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.key.size = unit(0.28, "lines"), legend.background = element_blank(),
        legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(2, 2, 1, 2))

att <- rbind(
  data.table(method = "squidpy\n(spatial)", class = hep_class(d$sq_source)),
  data.table(method = "LIANA\n(scRNA)",     class = hep_class(d$liana_source))
)[, .N, by = .(method, class)]
att[, method := factor(method, levels = c("LIANA\n(scRNA)", "squidpy\n(spatial)"))]
att[, class  := factor(class,  levels = c("Non-parenchymal", "Hepatocyte"))]
att_col <- c("Hepatocyte" = "#0D47A1", "Non-parenchymal" = "#C2185B")

p_bar <- ggplot(att, aes(method, N, fill = class)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = att_col, name = NULL,
                    guide = guide_legend(reverse = TRUE, keyheight = unit(0.3, "lines"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "LR pairs (source)") +
  theme_masld() +
  theme(plot.title = element_text(size = 6, hjust = 0, face = "plain"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(0.26, "lines"), legend.box.spacing = unit(1, "pt"),
        legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0),
        axis.text.x = element_text(size = 6),
        plot.margin = margin(2, 2, 1, 2))

# ==============================================================================
# COMPOSE — single row: scatter (concordance) + bar (source attribution)
# ==============================================================================
p <- p_sc + p_bar + plot_layout(widths = c(2.3, 1))

out <- file.path(FIG4_DIR, "panels", "figS4e.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, p, width = 5.2, height = 2.2, device = grDevices::cairo_pdf)
cat("[cccS4] saved:", out, "\n")

message(sprintf(paste0(
  "CAPTION (Fig S4e, spatial CCC consensus; demoted from main 4h): Spatial grounding of inferred cell-cell signaling (Visium spot-adjacency, ",
  "GSE192741). Of the ligand-receptor pairs inferred from dissociated single-cell RNA-seq (LIANA), ",
  "%d are independently re-detected in intact tissue by a spatial-proximity method (squidpy on Visium ",
  "spot adjacency). Left: each point is one of these %d shared pairs, positioned by its LIANA scRNA ",
  "disease score (x) and squidpy spatial Δ (y); teal = the %d/%d pairs (%.0f%%) whose disease direction ",
  "agrees between the two methods, grey = the rest. Right: the SOURCE cell type each method assigns to ",
  "the same %d pairs — squidpy attributes %d/%d to hepatocytes (%d/%d autocrine), LIANA only %d/%d to ",
  "hepatocytes and %d/%d to immune cells; a genuine methodological difference (spatial-proximity skews ",
  "hepatocyte-autocrine, dissociated-scRNA skews immune/stromal), not an error. The ~%.0f%% concordance ",
  "is reported honestly and is modest. PAIR-LEVEL recovery of inferred signaling, NOT cell-type-resolved ",
  "proof. Vu et al. Visium is deliberately excluded here: it has no healthy comparator and no per-spot ",
  "fibrosis staging, so it cannot support a disease-direction claim, and using GSE192741-healthy as its ",
  "reference would violate the established cross-cohort non-pooling rule. A companion single-cell-",
  "resolution exhibit (Govaere 2026 CosMx, IL32-hepatocyte -> macrophage) is reported separately as ",
  "Fig S4c."),
  n_pairs, n_pairs, n_conc, n_pairs, pct, n_pairs, sq_hep, n_pairs, sq_auto, n_pairs,
  liana_hep, n_pairs, liana_imm, n_pairs, pct))
