#!/usr/bin/env Rscript
# F5ef candidate Figure 5E: assay-native genetic-chromatin and protein-histology context for the
# complete 117-program frozen registry (R01-03, R11-03, G3-01).
#
# Changes from figures/main Fig 5E (style copied from scripts/figures/fig5e_multimodal_program_summary.R):
#   * universe: all 117 frozen v2 programs, not 22 programs selected on the superseded July disease fit.
#     Under the corrected fits no program set can be chosen by disease association (H59: 0 of 117),
#     so selection uses only assay-native results computed over the full 117 family.
#   * genetic column: member genes with >=1 adopted COLOC r2 primary pair (767 pairs; same release as
#     Fig 5D) whose SNP.PP.H4 posterior has >=0.5 mass in open chromatin of the program lineage
#     (01_genetic_context_117.py). Direct-trait (MASLD/liver fat) and liver-enzyme genes are stacked.
#   * histology column: rank-residual partial Spearman with participant bootstrap 95% CI; BH over the
#     testable programs of 117 (02_protein_histology_117.R).
#   * spatial column removed: its values used program-lineage adjustment, which differs from the
#     all-cell adjustment in Fig 5F, and an all-cell fit exists only for Stromal ECM and Ductular injury.
# Display rule: >=1 lineage-open colocalized gene OR robust histology association. Display filtering
# never re-tests; the source table carries all 117 rows.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
grDevices::pdf.options(useDingbats = FALSE)

P <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WS <- file.path(P, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef")
CAND <- file.path(P, "figures/candidates/review-r3-fixes-20261015T112148Z/F5ef")
dir.create(CAND, recursive = TRUE, showWarnings = FALSE)

FS <- 6
FAM <- "Helvetica"
INK <- "#222222"
NEG <- "#2166AC"
POS <- "#B2182B"
GOLD <- "#C79A3E"
GOLD_DARK <- "#7A5A1E"
NONROBUST <- "#C4C4C4"
# copied from ct_palette in scripts/figures/publication_theme.R
CT_COLS <- c(hepatocytes = "#0D47A1", cholangiocytes = "#1565C0", fibroblasts = "#F57F17",
             macrophages = "#C2185B", tcells = "#7B1FA2")
CT_ABBR <- c(hepatocytes = "Hep", cholangiocytes = "Chol", fibroblasts = "Fib",
             macrophages = "Mac", tcells = "T")
CT_ORDER <- names(CT_ABBR)

gen <- fread(file.path(WS, "out/genetic_context_117.tsv"))
hist <- fread(file.path(WS, "out/protein_histology_117.tsv"))
stopifnot(nrow(gen) == 117L, nrow(hist) == 117L, !anyDuplicated(gen$program_uid),
          setequal(gen$program_uid, hist$program_uid))
h <- hist[, .(program_uid, n_measured, retained_l1_weight, testable_scoring, testable, n, rho, ci_low,
              ci_high, pvalue, padj, sensitivity_sign_agree, batch_sign_agree, robust, correction_family)]
dt <- merge(gen, h, by = "program_uid", sort = FALSE)
# a gene open through both a direct-trait and a liver-enzyme pair is counted once, as direct
dt[, n_open_enzyme_only_genes := n_open_genes - n_open_direct_genes]
stopifnot(all(dt$n_open_enzyme_only_genes >= 0L))
dt[, display := n_open_genes > 0L | robust %in% TRUE]
dt[, display_rule := ">=1 lineage-open COLOC r2 gene OR robust protein-histology association"]
dt[, display_n := sum(display)]

sel <- dt[display == TRUE]
sel[, ct_rank := match(cell_type, CT_ORDER)]
setorder(sel, ct_rank, module)
sel[, y := .N:1L]
n <- nrow(sel)
sel[, short_name := sub(" \\([^()]+\\)$", "", program_name)]
sel[, row_label := paste0(module, "  ", short_name)]

out_tsv <- file.path(CAND, "fig5e_program_context_117_candidate.tsv")
fwrite(dt[, !c("open_pairs")], out_tsv, sep = "\t")

YMIN <- -1.45
YMAX <- n + 2.75
groups <- sel[, .(ymin = min(y) - 0.46, ymax = max(y) + 0.46, cy = mean(y),
                  colour = CT_COLS[.BY[[1]]], label = CT_ABBR[.BY[[1]]]), by = cell_type]
brk <- sel[, which(cell_type[-.N] != cell_type[-1L])]
breaks_y <- if (length(brk)) sel$y[brk] - 0.5 else numeric()

base_void <- theme_void(base_family = FAM, base_size = FS) +
  theme(text = element_text(family = FAM, size = FS, face = "plain", colour = INK),
        plot.margin = margin(1, 1, 1, 1), legend.position = "none")
separator_layers <- function(xmin, xmax) {
  if (!length(breaks_y)) return(list())
  list(geom_segment(data = data.table(y = breaks_y), aes(x = xmin, xend = xmax, y = y, yend = y),
                    inherit.aes = FALSE, linewidth = 0.22, linetype = "dashed", colour = "#9A9A9A"))
}

label_brackets <- rbind(
  groups[, .(x = 0.48, xend = 0.48, y = ymin + 0.12, yend = ymax - 0.12, colour)],
  groups[, .(x = 0.48, xend = 0.57, y = ymin + 0.12, yend = ymin + 0.12, colour)],
  groups[, .(x = 0.48, xend = 0.57, y = ymax - 0.12, yend = ymax - 0.12, colour)]
)
p_labels <- ggplot() +
  geom_segment(data = label_brackets, aes(x = x, xend = xend, y = y, yend = yend),
               colour = label_brackets$colour, linewidth = 0.35 / .pt, lineend = "round") +
  geom_text(data = groups, aes(x = 0.42, y = cy, label = label), hjust = 1, size = FS / .pt,
            family = FAM, colour = INK) +
  geom_text(data = sel, aes(x = 0.64, y = y, label = row_label), hjust = 0, size = FS / .pt,
            family = FAM, colour = INK) +
  separator_layers(0.03, 3.10) +
  coord_cartesian(xlim = c(0.02, 3.10), ylim = c(YMIN, YMAX), expand = FALSE) +
  base_void

# genetic-chromatin column: stacked direct-trait + liver-enzyme open colocalized genes
GEN_MAX <- max(3L, max(sel$n_open_genes))
x_ticks_gen <- data.table(x = 0:GEN_MAX)
p_genetic <- ggplot() +
  geom_rect(data = sel, aes(xmin = 0, xmax = GEN_MAX, ymin = y - 0.16, ymax = y + 0.16),
            fill = "#F1F1F1", colour = NA) +
  geom_rect(data = sel[n_open_direct_genes > 0L],
            aes(xmin = 0, xmax = n_open_direct_genes, ymin = y - 0.16, ymax = y + 0.16),
            fill = GOLD_DARK, colour = NA) +
  geom_rect(data = sel[n_open_enzyme_only_genes > 0L],
            aes(xmin = n_open_direct_genes, xmax = n_open_genes,
                ymin = y - 0.16, ymax = y + 0.16), fill = GOLD, colour = NA) +
  annotate("text", x = 0, y = c(n + 2.25, n + 1.55),
           label = c("Colocalized genes", "in open chromatin, n"), hjust = 0,
           size = FS / .pt, family = FAM, colour = INK) +
  annotate("rect", xmin = c(0, 0.48 * GEN_MAX), xmax = c(0.09 * GEN_MAX, 0.57 * GEN_MAX),
           ymin = n + 0.64, ymax = n + 0.86, fill = c(GOLD_DARK, GOLD), colour = NA) +
  annotate("text", x = c(0.12 * GEN_MAX, 0.60 * GEN_MAX), y = n + 0.75,
           label = c("direct", "enzyme"), hjust = 0, size = FS / .pt, family = FAM, colour = INK) +
  annotate("segment", x = 0, xend = GEN_MAX, y = -0.18, yend = -0.18, linewidth = 0.25, colour = INK) +
  geom_segment(data = x_ticks_gen, aes(x = x, xend = x, y = -0.18, yend = -0.32),
               linewidth = 0.25, colour = INK) +
  geom_text(data = x_ticks_gen, aes(x = x, y = -0.58, label = x), size = FS / .pt, family = FAM,
            colour = INK) +
  separator_layers(-0.08, GEN_MAX + 0.08) +
  coord_cartesian(xlim = c(-0.08, GEN_MAX + 0.08), ylim = c(YMIN, YMAX), expand = FALSE) +
  base_void

# protein-histology column: partial rho and participant-bootstrap 95% CI
hp <- copy(sel)
hp[, line_col := fifelse(robust %in% TRUE, fifelse(rho >= 0, POS, NEG), NONROBUST)]
HIST_LIM <- min(1, max(0.6, ceiling(max(abs(c(hp$ci_low, hp$ci_high)), na.rm = TRUE) * 10) / 10))
hist_ticks <- data.table(x = pretty(c(-HIST_LIM, HIST_LIM), n = 4L))[x >= -HIST_LIM & x <= HIST_LIM]
tst <- hp[testable == TRUE]
p_hist <- ggplot() +
  annotate("segment", x = 0, xend = 0, y = 0.48, yend = n + 0.48, linewidth = 0.2, colour = "#B8B8B8") +
  geom_segment(data = tst, aes(x = ci_low, xend = ci_high, y = y, yend = y),
               colour = tst$line_col, linewidth = 0.48, lineend = "round") +
  geom_segment(data = tst, aes(x = rho, xend = rho, y = y - 0.17, yend = y + 0.17),
               colour = tst$line_col, linewidth = 0.65, lineend = "round") +
  geom_text(data = hp[testable != TRUE], aes(x = 0.04, y = y, label = "not scored"), hjust = 0,
            size = FS / .pt, family = FAM, colour = INK) +
  annotate("text", x = 0, y = c(n + 2.25, n + 1.55),
           label = c("Protein vs histologic burden", "partial ρ (95% CI)"), hjust = 0.5,
           size = FS / .pt, family = FAM, colour = INK) +
  annotate("segment", x = c(-HIST_LIM, -0.86 * HIST_LIM, -0.02 * HIST_LIM),
           xend = c(-0.88 * HIST_LIM, -0.74 * HIST_LIM, 0.10 * HIST_LIM),
           y = n + 0.75, yend = n + 0.75, colour = c(POS, NEG, NONROBUST), linewidth = 0.48) +
  annotate("text", x = c(-0.70 * HIST_LIM, 0.14 * HIST_LIM), y = n + 0.75,
           label = c("robust", "not robust"), hjust = 0, size = FS / .pt, family = FAM, colour = INK) +
  annotate("segment", x = -HIST_LIM, xend = HIST_LIM, y = -0.18, yend = -0.18, linewidth = 0.25,
           colour = INK) +
  geom_segment(data = hist_ticks, aes(x = x, xend = x, y = -0.18, yend = -0.32), linewidth = 0.25,
               colour = INK) +
  geom_text(data = hist_ticks, aes(x = x, y = -0.58, label = sprintf("%g", x)), size = FS / .pt,
            family = FAM, colour = INK) +
  separator_layers(-HIST_LIM, HIST_LIM) +
  coord_cartesian(xlim = c(-HIST_LIM, HIST_LIM), ylim = c(YMIN, YMAX), expand = FALSE) +
  base_void

panel <- p_labels + p_genetic + p_hist + plot_layout(widths = c(1.85, 0.95, 1.30))
out_pdf <- file.path(CAND, "fig5e_program_context_117_candidate.pdf")
height_in <- 0.36 + 0.105 * (n + 4.2)
ggsave(out_pdf, panel, width = 3.7, height = height_in, units = "in", device = grDevices::cairo_pdf)

message("[F5ef 5E] ", out_pdf, "; displayed ", n, "/117; testable histology ", sum(dt$testable),
        "; q<0.05 ", sum(dt$padj < 0.05, na.rm = TRUE), "; robust ", sum(dt$robust),
        "; programs with >=1 open colocalized gene ", sum(dt$n_open_genes > 0L),
        "; with a direct-trait open gene ", sum(dt$n_open_direct_genes > 0L))
message("CAPTION (Fig. 5E candidate): Assay-native context for the complete 117-program frozen registry. ",
        "Rows show programs with at least one member gene in an adopted COLOC r2 primary pair whose ",
        "colocalization posterior lies mainly in open chromatin of the program's lineage, or a robust ",
        "protein-histology association; all 117 rows are in the source table. Histology estimates are ",
        "rank-residual partial Spearman correlations with participant-bootstrap 95% intervals (46 MASLD ",
        "participants), ", unique(dt$correction_family), ". Spatial organization is in Fig. 5F.")
writeLines(capture.output(sessionInfo()), file.path(WS, "out/sessionInfo_render_5e.txt"))
