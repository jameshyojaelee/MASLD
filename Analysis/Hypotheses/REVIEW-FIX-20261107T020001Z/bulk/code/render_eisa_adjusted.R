# KEY MESSAGE: Intronic and exonic stage effects align; adjusted ratio changes are highlighted.
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
ROOT <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/bulk'
B <- fread(file.path(ROOT, 'eisa_adjusted_v3/eisa_all_filtered_genes.tsv.gz'))
named <- c('GNMT', 'MAT1A', 'COL1A1', 'LUM', 'IGFBP7', 'TNFSF10')
B[, highlight := eisa_FDR < .05 & abs(dEmI) >= .5]
stopifnot(nrow(B) == 12010L, sum(B$eisa_FDR < .05) == 2711L, sum(B$highlight) == 114L)
fwrite(B[, .(gene, symbol, dE_log2 = dE, dI_log2 = dI, adjusted_log2_ratio_change = dEmI,
             ratio_BH_q = eisa_FDR, highlighted = highlight, class,
             ratio_model = '~ batch + inferred_sex + bin')],
       file.path(ROOT, 'figures/ssx_b_source.tsv'), sep = '\t')
FS <- 6; MAG <- '#C9265E'; BLU <- '#1565C0'; GRY <- '#9E9E9E'
g <- ggplot(B, aes(dE, dI)) +
  geom_abline(slope = 1, intercept = 0, linewidth = .25, colour = 'grey60', linetype = 2) +
  geom_point(data = B[highlight == FALSE], colour = GRY, size = .3, stroke = 0, alpha = .5) +
  geom_point(data = B[highlight == TRUE], colour = MAG, size = .6, stroke = 0) +
  geom_point(data = B[symbol %in% named], colour = BLU, size = 1, stroke = 0) +
  geom_text(data = B[symbol %in% named][, `:=`(hj = fifelse(symbol == 'COL1A1', 1.15, -.15),
               vj = fifelse(symbol == 'COL1A1', -.3, 1.1))],
            aes(label = symbol, hjust = hj, vjust = vj), colour = BLU, size = FS / .pt) +
  annotate('text', x = -Inf, y = Inf, hjust = -.05, vjust = 1.2, size = FS / .pt,
    label = sprintf('%s genes, Spearman %.2f\nmagenta: %d genes with adjusted exon/intron\nchange >= 1.41-fold, q < 0.05',
      format(nrow(B), big.mark = ','), cor(B$dE, B$dI, method = 'spearman'), sum(B$highlight))) +
  coord_equal() + labs(x = 'Exonic change, F3-F4 vs F0-F1 (log2)',
    y = 'Intronic change, F3-F4 vs F0-F1 (log2)', title = 'Intronic and exonic stage effects (493 donors)') +
  theme_classic(base_size = FS, base_family = 'Helvetica') +
  theme(text = element_text(size = FS, colour = 'black'), axis.text = element_text(size = FS, colour = 'black'),
    axis.title = element_text(size = FS), plot.title = element_text(size = FS, hjust = 0, margin = margin(b = 3)),
    axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3),
    plot.margin = margin(4, 6, 4, 4))
pdf.options(useDingbats = FALSE)
cairo_pdf(file.path(ROOT, 'figures/figs3_steatosite_exon_vs_intron.pdf'), width = 2.6, height = 2.7, family = 'Helvetica')
print(g); invisible(dev.off())
writeLines(capture.output(sessionInfo()), file.path(ROOT, 'figures/eisa_plot_sessionInfo.txt'))
