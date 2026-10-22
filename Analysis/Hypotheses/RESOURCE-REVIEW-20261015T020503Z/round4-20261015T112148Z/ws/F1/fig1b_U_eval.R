#!/usr/bin/env Rscript
# KEY MESSAGE: Among 3,769 genes evaluated by multi-signal COLOC, colocalized genes and
# canonical DEGs share 26, against 21.2 expected for genes of the same expression level.
# Copy of scripts/figures/fig1b_current_interface.R (layout kept), moved to the H57 U_eval
# universe per the H57 prespec null reading (RESULT_H57.md, run-20261015T113107Z).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(scales)})
grDevices::pdf.options(useDingbats = FALSE)
set.seed(20260923)  # no randomness used; set by project rule
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- file.path(root, "figures/candidates/review-r3-fixes-20261015T112148Z/F1")
setwd(root)
source("scripts/figures/publication_theme.R")
pdf_out <- file.path(out, "fig1b_complementary_maps_U_eval.pdf")
stopifnot(!file.exists(pdf_out))

g <- fread(file.path(out, "fig1b_U_eval_gene_table.tsv"))
N <- nrow(g); G <- sum(g$gen); D <- sum(g$deg); O <- sum(g$gen & g$deg)
stopifnot(N == 3769L, G == 426L, D == 220L, O == 26L)
# stratified expectation recomputed from the gene table (10 AveExpr bins, H57 "Bins")
E <- g[, .(e = sum(gen) * sum(deg) / .N), by = aveexpr_bin10][, sum(e)]
stopifnot(abs(E - 21.2127) < 1e-3)
never_evaluated <- 935L

view <- rbind(
  data.table(map = "Colocalized genes", segment = c("convergent", "genetic_only"),
             n = c(O, G - O)),
  data.table(map = "Canonical DEGs\n(BH q<0.05, |log2FC|>0.50)",
             segment = c("convergent", "disease_state_only"), n = c(O, D - O)))
view[, fraction := n / sum(n), by = map]
view[, segment := factor(segment, levels = c("convergent", "genetic_only", "disease_state_only"))]
setorder(view, map, segment)
view[, label := sprintf("%s\n%s", comma(n), percent(fraction, accuracy = .1))]
view[segment == "convergent", label := NA_character_]
lv <- c("Canonical DEGs\n(BH q<0.05, |log2FC|>0.50)", "Colocalized genes")
view[, map := factor(map, levels = lv)]
segs <- data.table(map = lv, x = c(E / D, E / G))
segs[, map := factor(map, levels = lv)]
shared <- merge(view[segment == "convergent", .(map, fraction, n)], segs[, .(map, xexp = x)], by = "map")
shared[, xlab := pmax(fraction, xexp)]
shared[, txt := sprintf("%s observed (%s)", comma(n), percent(fraction, accuracy = .1))]
note <- data.table(x = 0, y = 1 - .52,
  txt = sprintf("%s further canonical DEGs had no multi-signal COLOC result (not evaluated)", comma(never_evaluated)))
colors <- c(convergent = "#00695C", genetic_only = "#1565C0", disease_state_only = "#C9265E")
plot <- ggplot(view, aes(fraction, map, fill = segment)) +
  geom_col(width = .56, color = "white", linewidth = .35, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = label), position = position_stack(vjust = .5, reverse = TRUE),
            size = 6 / .pt, color = "white", lineheight = .9, na.rm = TRUE) +
  geom_segment(data = segs, inherit.aes = FALSE,
               aes(x = x, xend = x, y = as.numeric(map) - .40, yend = as.numeric(map) + .40),
               linetype = "22", linewidth = .4, color = "black") +
  geom_text(data = segs, inherit.aes = FALSE,
            aes(x = x, y = as.numeric(map) + .46,
                label = sprintf("expected at matched expression %s", percent(x, accuracy = .1))),
            size = 6 / .pt, color = "black", hjust = 0, nudge_x = .008) +
  geom_text(data = shared, inherit.aes = FALSE, aes(x = xlab, y = map, label = txt),
            size = 6 / .pt, color = "black", hjust = 0, nudge_x = .014) +
  geom_text(data = note, inherit.aes = FALSE, aes(x = x, y = y, label = txt),
            size = 6 / .pt, color = "black", hjust = 0) +
  scale_fill_manual(values = colors, name = NULL,
                    breaks = c("convergent", "disease_state_only", "genetic_only"),
                    labels = c("Both", "Disease-state only", "Colocalized only")) +
  scale_x_continuous(labels = percent, expand = expansion(mult = c(0, .02))) +
  scale_y_discrete(expand = expansion(add = .75)) +
  labs(x = "Fraction of each gene list", y = NULL) + theme_masld(base_size = 6) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
        legend.key.width = unit(.38, "cm"))
ggsave(pdf_out, plot, width = 4.8, height = 2.25, device = cairo_pdf, bg = "white")

summ <- data.table(universe = "U_eval (positive or tested_le_0.5 SuSiE state)", n = N,
  n_genetic = G, n_deg = D, n_both = O, n_genetic_only = G - O, n_deg_only = D - O,
  n_neither = N - G - D + O, n_deg_never_evaluated = never_evaluated,
  E_strat_10bin = E, O_over_E = O / E, marker_genetic_row = E / G, marker_deg_row = E / D,
  E_unstrat = G * D / N, unstrat_marker_genetic_row = (G * D / N) / G,
  unstrat_marker_deg_row = (G * D / N) / D,
  H57_MH_OR = 1.295, H57_MH_CI = "0.835-2.007", H57_exact_p = 0.2755,
  H57_label = "expression not fully controlled (S6)")
fwrite(summ, file.path(out, "fig1b_U_eval_summary.tsv"), sep = "\t")
fwrite(view, file.path(out, "fig1b_U_eval_display_source.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_fig1b.txt"))
message(sprintf("Caption: Among %d genes evaluated by multi-signal COLOC, %d colocalized and %d canonical DEGs share %d; %.2f expected for genes of the same expression level (10 AveExpr bins).", N, G, D, O, E))
