# KEY MESSAGES:
# 2D-ancestry: the direct-minus-enzyme coding/splice difference is large in EUR and near zero in non-EUR strata.
# 2D-subclass: the pooled direct excess is carried by the three PDFF GWAS; NAFLD/NASH alone has an interval that includes zero.
# companion: after liftover repair no noncoding posterior is unresolved; deposited-context fractions are similar between classes.
# 2A-support: most EUR enzyme genes sit outside any same-study fine-mapped window.
set.seed(20260923)
suppressPackageStartupMessages({library(ggplot2)})
grDevices::pdf.options(useDingbats = FALSE)
args <- commandArgs(TRUE); F <- args[1]
fs <- 6; col_d <- "#C9265E"; col_e <- "#1565C0"
th <- theme_classic(base_size = fs, base_family = "Helvetica") + theme(
  text = element_text(size = fs, colour = "black"), axis.text = element_text(size = fs, colour = "black"),
  axis.title = element_text(size = fs), legend.text = element_text(size = fs), legend.title = element_text(size = fs),
  plot.title = element_blank(), legend.key.size = unit(6, "pt"), axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25))
# Panel 1
a <- read.delim(file.path(F, "Fig2D_coding_by_ancestry_source.tsv"))
a$row <- factor(a$row, levels = rev(a$row))
a$lab <- sprintf("%.1f%% vs %.1f%%", 100*a$fraction_direct, 100*a$fraction_enzyme)
p1 <- ggplot(a, aes(y = row)) + geom_vline(xintercept = 0, colour = "#9E9E9E", linewidth = 0.25) +
  geom_segment(aes(x = 100*lo98_333, xend = 100*hi98_333, yend = row), colour = "grey60", linewidth = 0.3, na.rm = TRUE) +
  geom_segment(aes(x = 100*lo95, xend = 100*hi95, yend = row), colour = "black", linewidth = 0.5) +
  geom_point(aes(x = 100*difference), size = 1.2, shape = 21, fill = "black", colour = "black") +
  geom_text(aes(x = 52, label = lab), hjust = 0, size = fs/.pt) +
  scale_x_continuous("Direct minus enzyme, protein-altering + canonical-splice\nposterior probability (percentage points)", limits = c(-10, 75), breaks = c(-10, 0, 10, 20, 30, 40, 50)) +
  ylab(NULL) + th
ggsave(file.path(F, "Fig2D_coding_by_ancestry.pdf"), p1, width = 3.6, height = 1.5, device = cairo_pdf)
# Panel 2
s <- read.delim(file.path(F, "Fig2D_direct_subclass_source.tsv"))
s$row <- factor(s$row, levels = rev(s$row))
s$lab <- sprintf("%.1f%% vs %.1f%%; %d sets / %d clusters", 100*s$frac_direct, 100*s$frac_enzyme, s$n_cs_direct, s$n_clusters_direct)
p2 <- ggplot(s, aes(y = row)) + geom_vline(xintercept = 0, colour = "#9E9E9E", linewidth = 0.25) +
  geom_segment(aes(x = 100*diff_lo95, xend = 100*diff_hi95, yend = row), colour = "black", linewidth = 0.5) +
  geom_point(aes(x = 100*difference), size = 1.2, shape = 21, fill = "black") +
  geom_text(aes(x = 58, label = lab), hjust = 0, size = fs/.pt) +
  scale_x_continuous("Direct subset minus all liver enzymes,\nprotein-altering + canonical-splice posterior (percentage points)", limits = c(-5, 125), breaks = seq(0, 50, 10)) +
  ylab(NULL) + th
ggsave(file.path(F, "Fig2D_direct_subclass_split.pdf"), p2, width = 4.2, height = 1.15, device = cairo_pdf)
# Panel 3
x <- read.delim(file.path(F, "Fig2_companion_noncoding_context_repaired_source.tsv"))
x <- x[x$label == "credible-set-weighted repaired" & x$context != "unresolved", ]
lv <- c("promoter-proximal" = "Within 2 kb of TSS", "lineage-accessible" = "Open in a liver lineage", "ABC enhancer" = "Liver ABC enhancer", "resolved, no context" = "Resolved, no deposited context")
x$ctx <- factor(lv[x$context], levels = rev(lv))
xl <- rbind(data.frame(ctx = x$ctx, cls = "Direct MASLD / liver fat", v = 100*x$fraction_direct),
            data.frame(ctx = x$ctx, cls = "Liver enzymes", v = 100*x$fraction_enzyme))
p3 <- ggplot(xl, aes(x = v, y = ctx, colour = cls, shape = cls)) + geom_point(size = 1.3, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = sprintf("%.2f%%", v), x = v + 3), hjust = 0, size = fs/.pt, colour = "black", position = position_dodge(width = 0.5), show.legend = FALSE) +
  scale_colour_manual(NULL, values = c(col_d, col_e)) + scale_shape_manual(NULL, values = c(16, 17)) +
  scale_x_continuous("Share of noncoding posterior probability (%)", limits = c(0, 100), breaks = seq(0, 100, 25)) + ylab(NULL) + th +
  theme(legend.position = "top")
ggsave(file.path(F, "Fig2_companion_noncoding_context_repaired.pdf"), p3, width = 3.6, height = 1.7, device = cairo_pdf)
# Panel 4
c <- read.delim(file.path(F, "Fig2A_finemapped_support_by_cell_source.tsv"))
c <- c[c$multi_signal_genes_positive > 0, ]
c$cell <- paste(c$trait, c$ancestry)
c$cell <- factor(c$cell, levels = rev(c$cell[order(c$trait, c$ancestry)]))
cl <- rbind(data.frame(cell = c$cell, what = "Multi-signal genes with PP.H4 > 0.5 in this cell", n = c$multi_signal_genes_positive),
            data.frame(cell = c$cell, what = "of which lead inside a same-study fine-mapped window (+/-250 kb)", n = c$genes_in_same_study_window_pad250kb))
p4 <- ggplot(cl, aes(x = n, y = cell, shape = what)) + geom_point(size = 1.2, colour = "black") +
  geom_text(data = c, aes(x = 215, y = cell, label = finemapped_windows), inherit.aes = FALSE, size = fs/.pt, hjust = 1) +
  annotate("text", x = 215, y = nrow(c) + 0.9, label = "Fine-mapped windows", size = fs/.pt, hjust = 1, vjust = 0) +
  scale_shape_manual(NULL, values = c(1, 16)) + scale_x_continuous("Genes", limits = c(0, 215), breaks = seq(0, 200, 50)) +
  coord_cartesian(clip = "off") + ylab(NULL) + th + theme(legend.position = "bottom", legend.direction = "vertical", plot.margin = margin(10, 4, 2, 2))
ggsave(file.path(F, "Fig2A_finemapped_support_by_cell.pdf"), p4, width = 3.4, height = 2.6, device = cairo_pdf)
writeLines(capture.output(sessionInfo()), file.path(F, "sessionInfo_render_panels.txt"))
