# F2d repair panels (2026-10-15). Re-rendered into figures/candidates/.../F2d/repair/; first-round PDFs untouched.
# KEY MESSAGES:
# 2D-ancestry: the pooled 19.5-point coding/splice excess shrinks to 2.2 points (interval spans 0) in AMR, the one
#   ancestry with comparable fine-mapping coverage; no stratum receives a trait-class call.
# companion: after position re-annotation no noncoding posterior is unresolved; context fractions are similar
#   between classes when noncoding posterior is summed over all credible sets.
# 2A-support: 100 of 436 multi-signal genes have a lead inside a same-study fine-mapped window; the EUR enzyme cells
#   carry most genes and almost no fine-mapped windows.
set.seed(20260923)
suppressPackageStartupMessages({library(ggplot2)})
grDevices::pdf.options(useDingbats = FALSE)
args <- commandArgs(TRUE); F0 <- args[1]; F <- args[2]
fs <- 6; col_d <- "#C9265E"; col_e <- "#1565C0"
th <- theme_classic(base_size = fs, base_family = "Helvetica") + theme(
  text = element_text(size = fs, colour = "black"), axis.text = element_text(size = fs, colour = "black"),
  axis.title = element_text(size = fs), legend.text = element_text(size = fs), legend.title = element_text(size = fs),
  plot.title = element_blank(), legend.key.size = unit(6, "pt"), axis.line = element_line(linewidth = 0.25),
  axis.ticks = element_line(linewidth = 0.25))
fmt1 <- function(x) sprintf("%.1f", x)

# ---------------- Panel 1: coding/splice difference by ancestry stratum
a <- read.delim(file.path(F0, "Fig2D_coding_by_ancestry_source.tsv"))
stopifnot(nrow(a) == 5)
a$row_label <- c("Pooled, all ancestries (adopted interval)",
                 "EUR+AMR, fixed weights 16/26 and 10/26\nEUR enzyme arm: 8 of 767 lead windows fine-mapped",
                 "EUR only\nenzyme arm: 8 of 767 lead windows fine-mapped",
                 "AMR only",
                 "All non-EUR pooled")
a$status <- c("descriptive", "no call", "no call", "no call", "descriptive")
a$row_label <- factor(a$row_label, levels = rev(a$row_label))
a$lab <- sprintf("%s%% vs %s%%; %s (%s to %s); %s",
                 fmt1(100 * a$fraction_direct), fmt1(100 * a$fraction_enzyme),
                 fmt1(100 * a$difference), fmt1(100 * a$lo95), fmt1(100 * a$hi95), a$status)
p1 <- ggplot(a, aes(y = row_label)) + geom_vline(xintercept = 0, colour = "#9E9E9E", linewidth = 0.25) +
  geom_segment(aes(x = 100 * lo98_333, xend = 100 * hi98_333, yend = row_label), colour = "grey60",
               linewidth = 0.3, na.rm = TRUE) +
  geom_segment(aes(x = 100 * lo95, xend = 100 * hi95, yend = row_label), colour = "black", linewidth = 0.5) +
  geom_point(aes(x = 100 * difference), size = 1.2, shape = 21, fill = "black", colour = "black") +
  geom_text(aes(x = 54, label = lab), hjust = 0, size = fs / .pt) +
  scale_x_continuous("Direct minus enzyme, protein-altering + canonical-splice\nposterior probability (percentage points)",
                     breaks = c(-10, 0, 10, 20, 30, 40, 50)) +
  coord_cartesian(xlim = c(-10, 50), clip = "off") + ylab(NULL) + th + theme(plot.margin = margin(4, 175, 2, 2))
ggsave(file.path(F, "Fig2D_coding_by_ancestry_repair.pdf"), p1, width = 5.6, height = 1.8, device = cairo_pdf)
a_out <- a[, setdiff(names(a), "lab")]
a_out$row_label <- gsub("\n", "; ", as.character(a_out$row_label), fixed = TRUE)
a_out$printed_label <- gsub("\n", "; ", a$lab, fixed = TRUE)
write.table(a_out, file.path(F, "Fig2D_coding_by_ancestry_repair_source.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------- Panel 3: noncoding-context companion, estimator named
x <- read.delim(file.path(F0, "Fig2_companion_noncoding_context_repaired_source.tsv"))
x <- x[x$label == "credible-set-weighted repaired" & x$context != "unresolved", ]
stopifnot(nrow(x) == 4)
lv <- c("promoter-proximal" = "Within 2 kb of TSS", "lineage-accessible" = "Open in a liver lineage",
        "ABC enhancer" = "Liver ABC enhancer", "resolved, no context" = "Resolved, no deposited context")
x$ctx <- factor(lv[x$context], levels = rev(lv))
xl <- rbind(data.frame(ctx = x$ctx, cls = "Direct MASLD / liver fat", v = 100 * x$fraction_direct),
            data.frame(ctx = x$ctx, cls = "Liver enzymes", v = 100 * x$fraction_enzyme))
p3 <- ggplot(xl, aes(x = v, y = ctx, colour = cls, shape = cls)) +
  geom_point(size = 1.3, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = sprintf("%.2f%%", v), x = v + 3), hjust = 0, size = fs / .pt, colour = "black",
            position = position_dodge(width = 0.5), show.legend = FALSE) +
  scale_colour_manual(NULL, values = c(col_d, col_e)) + scale_shape_manual(NULL, values = c(16, 17)) +
  scale_x_continuous("Share of noncoding posterior probability (%),\nsummed over all credible sets of a class",
                     limits = c(0, 100), breaks = seq(0, 100, 25)) + ylab(NULL) + th +
  theme(legend.position = "top")
ggsave(file.path(F, "Fig2_companion_noncoding_context_repair.pdf"), p3, width = 3.6, height = 1.8, device = cairo_pdf)
xs <- x; xs$estimator <- "sum of context posterior over sum of noncoding posterior, all gated credible sets of the class"
write.table(xs, file.path(F, "Fig2_companion_noncoding_context_repair_source.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

# ---------------- Panel 4: Figure 2A support, all 20 cells and the totals
c4 <- read.delim(file.path(F0, "Fig2A_finemapped_support_by_cell_source.tsv"))
stopifnot(nrow(c4) == 20, sum(c4$finemapped_windows) == 1665)
c4$cell <- paste(c4$trait, c4$ancestry)
ord <- c4$cell[order(c4$trait, c4$ancestry)]
tot_lab <- "All cells, unique genes"
c4$cell <- factor(c4$cell, levels = rev(c(ord, tot_lab)))
cl <- rbind(data.frame(cell = c4$cell, what = "Multi-signal genes with PP.H4 > 0.5 in this cell",
                       n = c4$multi_signal_genes_positive),
            data.frame(cell = c4$cell, what = "of which lead inside a same-study fine-mapped window (+/-250 kb)",
                       n = c4$genes_in_same_study_window_pad250kb))
num <- data.frame(cell = factor(c(as.character(c4$cell), tot_lab), levels = levels(c4$cell)),
                  g = c(c4$multi_signal_genes_positive, 436), w = c(c4$genes_in_same_study_window_pad250kb, 100),
                  f = c(c4$finemapped_windows, 1665))
p4 <- ggplot(cl, aes(x = n, y = cell, shape = what)) + geom_point(size = 1.2, colour = "black") +
  geom_text(data = num, aes(x = 245, y = cell, label = g), inherit.aes = FALSE, size = fs / .pt, hjust = 1) +
  geom_text(data = num, aes(x = 290, y = cell, label = w), inherit.aes = FALSE, size = fs / .pt, hjust = 1) +
  geom_text(data = num, aes(x = 335, y = cell, label = f), inherit.aes = FALSE, size = fs / .pt, hjust = 1) +
  annotate("text", x = c(245, 290, 335), y = nrow(num) + 0.9, label = c("Genes", "In window", "Windows"),
           size = fs / .pt, hjust = 1, vjust = 0) +
  scale_shape_manual(NULL, values = c(1, 16)) +
  scale_x_continuous("Genes", limits = c(0, 335), breaks = seq(0, 200, 50)) +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(clip = "off") + ylab(NULL) + th +
  theme(legend.position = "bottom", legend.direction = "vertical", plot.margin = margin(10, 4, 2, 2))
ggsave(file.path(F, "Fig2A_finemapped_support_by_cell_repair.pdf"), p4, width = 4.0, height = 3.0, device = cairo_pdf)
write.table(num, file.path(F, "Fig2A_finemapped_support_by_cell_repair_source.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(F, "sessionInfo_render_repair.txt"))
