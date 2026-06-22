#!/usr/bin/env Rscript
# fig2_coding_noncoding.R  (2026-06-12)  — Fig 2C
# Where does the colocalizing causal variant fall? Coding vs non-coding split.
#
# Answers the draft's Fig 2C placeholder "[How many noncoding vs coding variants?]".
# Stacked by the VEP CONSEQUENCE of the COLOC top-H4 SNP (NOT a TSS-distance cut),
# for the two coloc methods. Headline: only 31/751 (4.1%) colocalizing genes are
# coding-led; ~96% act through non-coding regulatory sequence -> motivates the
# expression-mediated effector model (and the coding class deferred to later figs).
#
# Data: RNA-seq/results/coloc_variant_classes/{panel_F1_data.csv, panel_F2_hypergeom.csv}
# Out : figures/main/fig2_genetics/panels/coding_noncoding_split.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")
DAT <- file.path(BASE, "RNA-seq/results/coloc_variant_classes")

f1 <- fread(file.path(DAT, "panel_F1_data.csv"))
hg <- fread(file.path(DAT, "panel_F2_hypergeom.csv"))

# consequence order — ggplot stacks the FIRST level at the TOP, so put the coding
# classes first to cap each bar (makes the tiny coding fraction read at a glance)
lev <- c("coding","spliceSite","fiveUTR","threeUTR","promoter","intron","intergenic")
labs <- c(intergenic="intergenic", intron="intron", promoter="promoter",
          threeUTR="3' UTR", fiveUTR="5' UTR", spliceSite="splice site", coding="coding")
f1[, fine_class := factor(fine_class, levels = lev)]
f1[, method := factor(method, levels = c("ABF-COLOC","SuSiE-COLOC"))]

cls_cols <- c(intergenic="#B0BEC5", intron="#64B5F6", promoter="#1565C0",
              threeUTR="#4DB6AC", fiveUTR="#80CBC4",
              spliceSite="#880E4F", coding="#C9265E")

tot  <- f1[, .(n = sum(n_genes)), by = method]
codg <- f1[fine_class %in% c("coding","spliceSite"), .(coding = sum(n_genes)), by = method]
ann  <- merge(tot, codg, by = "method")
ann[, pct := round(100 * coding / n, 1)]
all_n <- hg$n_all_coloc[1]; all_cod <- hg$n_all_coloc_coding[1]

p <- ggplot(f1, aes(x = method, y = n_genes, fill = fine_class)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.15) +
  geom_text(data = tot, aes(x = method, y = n + 14, label = n),
            inherit.aes = FALSE, size = 2.6, fontface = "bold") +
  # call out the (tiny) coding-led fraction at the top of each bar
  geom_text(data = ann, aes(x = method, y = n - coding/2, label = paste0(coding, " coding")),
            inherit.aes = FALSE, size = 1.9, color = "white", fontface = "bold") +
  scale_fill_manual(values = cls_cols, labels = labs, name = "Consequence of\nCOLOC lead SNP",
                    breaks = lev) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Colocalizing genes (PP.H4 > 0.5)",
       title = "Colocalizing variants are overwhelmingly non-coding",
       subtitle = sprintf("Only %d/%d (%.1f%%) colocalizing genes are coding-led — ~96%% act through non-coding sequence",
                          all_cod, all_n, 100*all_cod/all_n)) +
  theme_masld() + theme_pub() +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = PUB_SUBTITLE, color = "grey30"),
        axis.text.x = element_text(size = PUB_AXIS_TITLE, face = "bold"))

save_fig(p, file.path(PANEL_DIR, "coding_noncoding_split.pdf"),
         width = fig_half_width + 0.6, height = 3.4)

out <- merge(f1[, .(method, consequence = as.character(fine_class), coarse_class, n_genes)],
             ann[, .(method, method_total = n, coding_led = coding, coding_pct = pct)],
             by = "method")
fwrite(out, file.path(PANEL_DIR, "coding_noncoding_source.csv"))
cat(sprintf("[fig2C] wrote coding_noncoding_split.pdf  (all-coloc coding-led = %d/%d = %.1f%%)\n",
            all_cod, all_n, 100*all_cod/all_n))
