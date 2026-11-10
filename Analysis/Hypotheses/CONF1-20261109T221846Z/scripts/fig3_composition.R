# CONF1 candidate panel (not adopted). Individual PDF, 6 pt, no titles or stats in the panel.
# KEY MESSAGE: discovery F3-F4 genes change in the same direction and more strongly in SteatoSITE,
#              and adjusting for deconvolved cell fractions removes most of that change.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE); OUT <- args[1]
if (is.na(OUT) || dir.exists(OUT)) stop("output directory missing or already exists")
dir.create(OUT); pdf.options(useDingbats = FALSE); set.seed(20261109)
G <- fread("out_c3_partition/gene_partition.tsv.gz")
G <- G[is.finite(D0_F3F4_log2FC) & is.finite(bulk_log2fc) & is.finite(adj_log2fc)]
FS <- 6
th <- theme_classic(base_size = FS, base_family = "Helvetica") +
  theme(axis.text = element_text(size = FS, colour = "black"), axis.title = element_text(size = FS),
        strip.text = element_text(size = FS), strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25))
L <- rbind(G[, .(symbol, x = D0_F3F4_log2FC, y = bulk_log2fc, fit = "Unadjusted")],
           G[, .(symbol, x = D0_F3F4_log2FC, y = adj_log2fc, fit = "Cell-fraction adjusted")])
L[, fit := factor(fit, levels = c("Unadjusted", "Cell-fraction adjusted"))]
sl <- L[, .(n = .N, slope_lm = coef(lm(y ~ x))[2], slope_origin = coef(lm(y ~ 0 + x))[2],
            same_sign = mean(sign(x) == sign(y))), by = fit]
lim <- quantile(abs(c(L$x, L$y)), 0.999)
p <- ggplot(L, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_point(colour = "#9E9E9E", size = 0.25, alpha = 0.35, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#C9265E", linewidth = 0.4) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) + facet_wrap(~fit, nrow = 1) +
  labs(x = "Discovery log2FC, F3-F4 vs F0-F1", y = "SteatoSITE log2FC, F3-F4 vs F0-F1") + th
ggsave(file.path(OUT, "fig3cand_C_steatosite_composition.pdf"), p, width = 3.2, height = 1.7, units = "in", device = cairo_pdf)
fwrite(sl, file.path(OUT, "fig3cand_C_values.tsv"), sep = "\t"); print(sl)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
