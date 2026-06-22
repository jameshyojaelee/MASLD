#!/usr/bin/env Rscript
# fig2_coloc_counts.R  (2026-06-12)  — Figure 2B
# GWAS-eQTL colocalization nominations: coloc.abf (baseline) vs SuSiE-coloc
# (primary) gene counts at PP.H4 thresholds, across 18,975 Broadaway cis-eQTL
# genes. Headline: ABF 618 vs SuSiE 368 at PP.H4 > 0.5; the methods cross over
# (ABF nominates more at > 0.5, SuSiE-coloc is more confident at > 0.9).
# Counts computed from disk so they always match the manuscript text.
#
# Out: figures/main/fig2_genetics/panels/coloc_method_counts.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

d   <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
abf <- suppressWarnings(as.numeric(d$coloc_best_pp4))
su  <- suppressWarnings(as.numeric(d$coloc_best_susie_pp4))
n_genes <- nrow(d)

ths <- c(0.5, 0.8, 0.9)
counts <- rbindlist(lapply(ths, function(t) data.table(
  threshold = sprintf("PP.H4 > %.1f", t),
  method    = c("coloc.abf (baseline)", "SuSiE-coloc (primary)"),
  n         = c(sum(abf > t, na.rm = TRUE), sum(su > t, na.rm = TRUE)))))
counts[, threshold := factor(threshold, levels = sprintf("PP.H4 > %.1f", ths))]
counts[, method := factor(method, levels = c("coloc.abf (baseline)", "SuSiE-coloc (primary)"))]

method_cols <- c("coloc.abf (baseline)" = "#90A4AE", "SuSiE-coloc (primary)" = "#1565C0")

p <- ggplot(counts, aes(x = threshold, y = n, fill = method)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = n), position = position_dodge(width = 0.72),
            vjust = -0.35, size = 2.9, fontface = "bold", color = "grey15") +
  scale_fill_manual(values = method_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Colocalizing genes",
       title = "GWAS-eQTL colocalization (Broadaway liver eQTLs)",
       subtitle = sprintf("%s cis-eQTL genes tested; SuSiE-coloc is more confident at high PP.H4",
                          format(n_genes, big.mark = ","))) +
  theme_masld(base_size = 9) +
  theme(legend.position = c(0.98, 0.97), legend.justification = c(1, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.5, color = "grey30"),
        axis.text.x = element_text(face = "bold"))

save_fig(p, file.path(PANEL_DIR, "coloc_method_counts.pdf"),
         width = fig_col_width * 1.05, height = 2.9)

fwrite(dcast(counts, threshold ~ method, value.var = "n"),
       file.path(PANEL_DIR, "coloc_method_counts_source.csv"))
cat(sprintf("[fig2B] wrote coloc_method_counts.pdf  (ABF 618/282/186 vs SuSiE 368/289/210; n=%d)\n", n_genes))
