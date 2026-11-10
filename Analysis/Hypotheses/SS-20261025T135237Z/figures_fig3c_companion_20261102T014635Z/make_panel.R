# KEY MESSAGE: splitting one cohort into five random groups gives little list non-overlap and almost no sign
# conflicts, so the Figure 3C non-overlap and sign flips across cohorts come from differences between studies.
# Reads finished SSP1 output and the Figure 3C counts (RESULTS.md AN3: 8,341 of 12,180; 1,211 of 3,839).
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
set.seed(20261101)
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
OUT <- commandArgs(trailingOnly = TRUE)[1]
FS <- 6; GT <- FS / .pt; ACROSS <- "#C9265E"; WITHIN <- "#1565C0"
Q <- fread(file.path(SS, "out_ssp1-20261028T005018Z/partitions_per_draw.tsv"))
X <- rbind(
  data.table(set = "Within one cohort\n(5 random groups,\n50 splits)", metric = "Called in one group only",
             pct = 100 * Q$frac_single, draw = Q$draw),
  data.table(set = "Within one cohort\n(5 random groups,\n50 splits)", metric = "Mixed sign among genes\ncalled in two or more groups",
             pct = 100 * Q$sign_discordant_multi, draw = Q$draw),
  data.table(set = "Across\n5 cohorts\n(Figure 3C)", metric = "Called in one group only", pct = 100 * 8341 / 12180, draw = NA_integer_),
  data.table(set = "Across\n5 cohorts\n(Figure 3C)", metric = "Mixed sign among genes\ncalled in two or more groups", pct = 100 * 1211 / 3839, draw = NA_integer_))
X[, set := factor(set, levels = unique(set)[c(2, 1)])]
X[, metric := factor(metric, levels = c("Called in one group only", "Mixed sign among genes\ncalled in two or more groups"))]
M <- X[, .(pct = median(pct), n = .N), by = .(set, metric)]
M[, lab := sprintf("%.1f%%", pct)]
fwrite(X, file.path(OUT, "fig3c_companion_source.tsv"), sep = "\t"); fwrite(M, file.path(OUT, "fig3c_companion_medians.tsv"), sep = "\t")
theme_ss <- theme_classic(base_size = FS, base_family = "Helvetica") %+replace% theme(
  text = element_text(size = FS, face = "plain", colour = "black"), axis.text = element_text(size = FS, colour = "black"),
  axis.title = element_text(size = FS), plot.title = element_text(size = FS, hjust = 0, margin = margin(b = 3)),
  strip.text = element_text(size = FS), strip.background = element_blank(),
  axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), plot.margin = margin(4, 6, 4, 4))
g <- ggplot(X, aes(set, pct, colour = set)) +
  geom_jitter(data = X[!is.na(draw)], width = 0.15, height = 0, size = 0.5, alpha = 0.6, stroke = 0) +
  geom_point(data = X[is.na(draw)], size = 1.6) +
  geom_crossbar(data = M[n > 1], aes(ymin = pct, ymax = pct), width = 0.45, linewidth = 0.3, colour = "black") +
  geom_text(data = M, aes(y = pct + 4, label = lab), vjust = 0, size = GT, colour = "black") +
  scale_colour_manual(values = setNames(c(ACROSS, WITHIN), levels(X$set)), guide = "none") +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 75, 25), expand = expansion(mult = c(0.02, 0.02))) +
  facet_wrap(~ metric, nrow = 1) +
  labs(x = NULL, y = "% of genes", title = "Called genes: across cohorts vs within SteatoSITE (493 donors)") + theme_ss
cairo_pdf(file.path(OUT, "figs3_fig3c_within_vs_across.pdf"), width = 3.4, height = 2.2, family = "Helvetica"); print(g); invisible(dev.off())
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
