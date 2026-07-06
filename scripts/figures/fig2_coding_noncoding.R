#!/usr/bin/env Rscript
# fig2_coding_noncoding.R  (2026-06-12)  — Fig 2C
# Where does the colocalizing causal variant fall? Coding vs non-coding split.
#
# Answers the draft's Fig 2C placeholder "[How many noncoding vs coding variants?]".
# Single horizontal stacked bar over the de-duplicated SuSiE-OR-ABF colocalization set
# (best PP.H4 > 0.5 per gene), segmented by the VEP CONSEQUENCE of that gene's strongest
# COLOC lead SNP (NOT a TSS-distance cut). Headline: only 59/1,527 (3.9%) colocalizing genes
# are coding-led; ~96% act through non-coding regulatory sequence -> motivates the
# expression-mediated effector model (and the coding class deferred to later figs).
#
# 2026-06-25: switched from SuSiE-only (736) to the de-duped SuSiE-OR-ABF union (1,527) so the
#             denominator MATCHES Fig2G (ancestry specificity, also SuSiE/ABF). Per-gene
#             consequence = the fine_class of the gene's max-PP.H4 colocalization.
# 2026-07-04: portfolio rebuilt with MVP (23->50 GWAS); union grew 751->1,527, coding 26->59.
#             All counts recomputed from disk, so the plot tracks the rebuild automatically.
# Data: RNA-seq/results/coloc_variant_classes/coloc_variant_annotation.csv (per-gene/GWAS,
#       carries pp4_susie/pp4_abf/pp4_best + fine_class/coarse_class).
# Out : figures/main/fig2_genetics/panels/Fig2C_coding_noncoding_split.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")
DAT <- file.path(BASE, "RNA-seq/results/coloc_variant_classes")

# de-duped SuSiE U ABF: one row per gene = its strongest colocalization (max best PP.H4)
av <- fread(file.path(DAT, "coloc_variant_annotation.csv"))[!is.na(pp4_best)]
setorder(av, gene_symbol, -pp4_best)
best <- av[, .SD[1], by = gene_symbol][pp4_best > 0.5]
# Merge the single spliceSite gene (1 gene) into coding to simplify the legend and bar
best[fine_class == "spliceSite", fine_class := "coding"]
f1 <- best[, .(n_genes = .N), by = fine_class]
f1[, `:=`(method = "SuSiE or ABF",                 # single bar over the union set
          coarse_class = fifelse(fine_class == "coding", "coding", "non-coding"))]

# consequence order — coding first so that coding lands at the very top of the vertical stack
lev <- c("coding","fiveUTR","threeUTR","promoter","intron","intergenic")
labs <- c(intergenic="intergenic", intron="intron", promoter="promoter",
          threeUTR="3' UTR", fiveUTR="5' UTR", coding="coding")
f1[, fine_class := factor(fine_class, levels = lev)]
f1[, method := factor(method, levels = "SuSiE or ABF")]
present_lev <- lev[lev %in% as.character(f1$fine_class)]   # legend shows only classes with data

cls_cols <- c(intergenic="#B0BEC5", intron="#64B5F6", promoter="#1565C0",
              threeUTR="#4DB6AC", fiveUTR="#80CBC4",
              spliceSite="#880E4F", coding="#C9265E")

tot  <- f1[, .(n = sum(n_genes)), by = method]
codg <- f1[fine_class %in% c("coding","spliceSite"), .(coding = sum(n_genes)), by = method]
ann  <- merge(tot, codg, by = "method")
ann[, pct := round(100 * coding / n, 1)]

ntot <- ann$n; ncod <- ann$coding; nnon <- ntot - ncod

p <- ggplot(f1, aes(x = method, y = n_genes, fill = fine_class)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.18) +
  # counts inside the wide segments (white); thin slivers are named by the coarse tags below
  geom_text(aes(label = ifelse(n_genes >= 30, n_genes, "")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 2.7) +
  # total above the vertical bar (black)
  geom_text(data = tot, aes(x = method, y = n, label = paste0(n, " genes")),
            inherit.aes = FALSE, vjust = -0.4, hjust = 0.5, size = 2.7, fontface = "bold") +
  # coarse coding-vs-non-coding headline, next to the bar (black)
  # non-coding: a bracket spanning the regulatory classes, label centred to the right (bracket opens leftward)
  annotate("segment", x = 1.34, xend = 1.34, y = 5, yend = nnon - 5, linewidth = 0.4, colour = "black") +
  annotate("segment", x = 1.34, xend = 1.30, y = 5, yend = 5, linewidth = 0.4, colour = "black") +
  annotate("segment", x = 1.34, xend = 1.30, y = nnon - 5, yend = nnon - 5, linewidth = 0.4, colour = "black") +
  annotate("text", x = 1.38, y = nnon/2, hjust = 0, vjust = 0.5, size = 2.7, fontface = "bold", colour = "black",
           label = sprintf("non-coding  %d (%.0f%%)", nnon, 100 - ann$pct)) +
  # coding: a leader line to the (tiny) coding sliver, label to the right (line is black)
  annotate("segment", x = 1.26, xend = 1.36, y = nnon + ncod/2, yend = nnon + ncod/2,
           linewidth = 0.4, colour = "black") +
  annotate("text", x = 1.38, y = nnon + ncod/2, hjust = 0, vjust = 0.5, size = 2.5,
           fontface = "bold", colour = "black",
           label = sprintf("coding  %d (%.0f%%)", ncod, ann$pct)) +
  scale_fill_manual(values = cls_cols, labels = unname(labs[present_lev]),
                    name = NULL, breaks = present_lev) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  scale_x_discrete(expand = expansion(add = c(0.45, 0.75))) +
  labs(x = NULL, y = "Colocalizing genes (best PP.H4 > 0.5)",
       title = "Colocalizing variant class") +
  coord_cartesian(clip = "off") +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom", legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 7),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
  guides(fill = guide_legend(ncol = 2))

save_fig(p, file.path(PANEL_DIR, "Fig2C_coding_noncoding_split.pdf"),
         width = 2.5, height = 3.8)

out <- merge(f1[, .(method, consequence = as.character(fine_class), coarse_class, n_genes)],
             ann[, .(method, method_total = n, coding_led = coding, coding_pct = pct)],
             by = "method")
fwrite(out, file.path(PANEL_DIR, "Fig2C_coding_noncoding_source.csv"))
cat(sprintf("[fig2C] wrote Fig2C_coding_noncoding_split.pdf  (SuSiE-or-ABF coding-led = %d/%d = %.1f%%)\n",
            ann$coding, ann$n, ann$pct))
message(sprintf("CAPTION (Fig2C): Consequence class of each colocalizing gene's strongest COLOC lead SNP (best PP.H4 > 0.5, SuSiE or ABF; de-duplicated; same %d-gene set as Fig2G). Only %d (%.1f%%) are coding-led; %.0f%% act through non-coding regulatory sequence.",
                ann$n, ann$coding, ann$pct, 100 - ann$pct))
