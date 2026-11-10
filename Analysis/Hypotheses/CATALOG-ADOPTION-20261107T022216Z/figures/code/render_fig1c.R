#!/usr/bin/env Rscript
# KEY MESSAGE: Of the 436 genes colocalized with liver eQTL (SuSiE PP.H4 > 0.5), 424 are
# supported in ALT/AST/GGT GWAS and 22 in NAFLD/NASH/PDFF GWAS; 368 of the 424 were never
# SuSiE-tested in a NAFLD/NASH/PDFF GWAS, so "enzyme only" is established for 46 genes.
#
# Copied from scripts/figures/fig1c_corrected_genetics.R (producer named in
# figures/main/fig1_atlas_overview/fig1c_input_manifest.tsv) and redrawn to the 1C redraw spec
# of Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/PANEL_SPECS.md, with recommended decision 1
# (REPORT.md): "not SuSiE-tested in the other trait class" is drawn as its own outlined state.
# Changes: three bars (enzyme only / both / direct only) become one stacked bar per trait class,
# split by the gene's SuSiE state in the other class (supported, tested at PP.H4 <= 0.5, not
# SuSiE-tested); tested denominators, study counts, the 436 total and the 362-gene common
# denominator are printed; fonts embedded via cairo_pdf; fills no longer reuse 1B or Liang
# semantic colours.
#
# Values: read from the deep-review cross-tab (recompute_1c_crosstab.json, built from
# GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv, 35 main
# tier-1/2 studies) and checked against recompute_1c.json and the adopted 414/10/12 partition
# (figures/main/fig1_atlas_overview/fig1c_trait_scope_source.tsv).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(jsonlite); library(scales)
})
grDevices::pdf.options(useDingbats = FALSE)
set.seed(20261017)  # no randomness used; set by project rule
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
setwd(root)
source("scripts/figures/publication_theme.R")
out <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/figures"
pdf_out <- file.path(out, "panels/fig1c_genetic_trait_scope.pdf")
src_out <- file.path(out, "source_tables/fig1c.tsv")
stopifnot(!file.exists(pdf_out), !file.exists(src_out))

deep <- file.path(root, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1")
ct <- fromJSON(file.path(deep, "recompute_1c_crosstab.json"))
rc <- fromJSON(file.path(deep, "recompute_1c.json"))
ut <- fromJSON(file.path(deep, "recompute_1c_untested.json"))
adopted <- fread("figures/main/fig1_atlas_overview/fig1c_trait_scope_source.tsv")
# pandas DataFrame.to_dict(): ct$crosstab[[direct_state]][[enzyme_state]]
# (outer key = direct-trait SuSiE state, inner key = enzyme SuSiE state; checked by margins)
x <- ct$crosstab
stopifnot(x$All$positive == 424L, x$positive$All == 22L, x$All$All == 18976L)
enz_pos <- c(supported = x$positive$positive, tested = x$tested_le_0.5$positive, none = x$none$positive)
dir_pos <- c(supported = x$positive$positive, tested = x$positive$tested_le_0.5, none = x$positive$none)
stopifnot(all(unname(enz_pos) == c(10, 46, 368)), all(unname(dir_pos) == c(10, 11, 1)),
          sum(enz_pos) == rc$liver_enzyme_positive_genes, sum(dir_pos) == rc$direct_MASLD_positive_genes,
          sum(enz_pos) + sum(dir_pos) - enz_pos[["supported"]] == rc$gene_level_positive,
          enz_pos[["tested"]] + enz_pos[["none"]] == adopted[genetic_trait_scope == "enzyme", n_genes],
          dir_pos[["tested"]] + dir_pos[["none"]] == adopted[genetic_trait_scope == "direct_disease", n_genes],
          enz_pos[["supported"]] == adopted[genetic_trait_scope == "both", n_genes],
          enz_pos[["none"]] == ut$enzyme_only$untested_in_other_class,
          dir_pos[["none"]] == ut$direct_only$untested_in_other_class)
n_tot <- rc$gene_level_positive
n_enz_studies <- rc$main_by_class$liver_enzyme; n_dir_studies <- rc$main_by_class$direct_MASLD
tested_enz <- rc$liver_enzyme_tested_genes; tested_dir <- rc$direct_MASLD_tested_genes
both_tested <- ct$n_tested_both_classes
both_enz <- ct$both_tested_enzyme_pos; both_dir <- ct$both_tested_direct_pos
stopifnot(n_tot == 436L, n_enz_studies == 20L, n_dir_studies == 15L, tested_enz == 3888L, tested_dir == 392L,
          both_tested == 362L, both_enz == 56L, both_dir == 21L,
          both_enz == enz_pos[["supported"]] + enz_pos[["tested"]],
          both_dir == dir_pos[["supported"]] + dir_pos[["tested"]])

state_lab <- c(supported = "Also supported in the other trait class",
               tested = "Finite result in the other class, PP.H4 ≤ 0.5",
               none = "No finite result in the other class")
row_lab <- c(enzyme = sprintf("Colocalized in ALT/AST/GGT\nGWAS (%d studies)", n_enz_studies),
             direct = sprintf("Colocalized in NAFLD/NASH/PDFF\nGWAS (%d studies)", n_dir_studies))
mk <- function(cls, v, y, tested, nst, other) {
  data.table(trait_class = cls, other_class = other, y = y, state = names(v), n = unname(v),
             x_start = cumsum(c(0, head(unname(v), -1))), x_end = cumsum(unname(v)),
             class_total = sum(v), class_tested_genes = tested, class_main_studies = nst)
}
seg <- rbind(mk("liver_enzyme (ALT/AST/GGT)", enz_pos, 2, tested_enz, n_enz_studies, "direct (NAFLD/NASH/PDFF)"),
             mk("direct (NAFLD/NASH/PDFF)", dir_pos, 1, tested_dir, n_dir_studies, "liver_enzyme (ALT/AST/GGT)"))
seg[, state_label := state_lab[state]]
seg[, crosstab_cell := fifelse(trait_class %like% "^liver",
                               sprintf("enzyme=positive, direct=%s", c(supported = "positive", tested = "tested_le_0.5", none = "none")[state]),
                               sprintf("direct=positive, enzyme=%s", c(supported = "positive", tested = "tested_le_0.5", none = "none")[state]))]
seg[, source := "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1/recompute_1c_crosstab.json"]
fills <- c(supported = "#6470d0", tested = "#aca5ed", none = "white")
seg[, fill := fills[state]]
seg[, outline := fifelse(state == "none", "black", "none")]
txt <- seg[, .(y = y[1] - 0.40, label = sprintf("%d = %d + %d + %d, of %s finite SuSiE results",
                                                class_total[1], n[1], n[2], n[3], comma(class_tested_genes[1]))),
           by = trait_class]
seg[, row_end_label := txt$label[match(trait_class, txt$trait_class)]]
fwrite(seg, src_out, sep = "\t")
extra <- data.table(trait_class = "both classes tested", quantity = c("genes SuSiE-tested in both classes",
  "enzyme-supported among them", "direct-supported among them", "supported in both", "total colocalized genes"),
  n = c(both_tested, both_enz, both_dir, enz_pos[["supported"]], n_tot),
  source = c(rep("recompute_1c_crosstab.json", 4), "recompute_1c.json gene_level_positive"))
fwrite(extra, file.path(out, "source_tables/fig1c_common_denominator.tsv"), sep = "\t")

seg[, state := factor(state, levels = names(state_lab))]
plot <- ggplot(seg) +
  geom_rect(aes(xmin = x_start, xmax = x_end, ymin = y - 0.27, ymax = y + 0.27, fill = state, colour = state),
            linewidth = 0.3) +
  scale_fill_manual(values = fills, labels = state_lab, name = NULL) +
  scale_colour_manual(values = c(supported = NA, tested = NA, none = "black"), labels = state_lab, name = NULL) +
  geom_text(data = txt, aes(x = 4, y = y, label = label), hjust = 0, vjust = 1, size = 6 / .pt, colour = "black") +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 400, 100), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(2, 1), labels = unname(row_lab[c("enzyme", "direct")]),
                     limits = c(0.45, 2.35), expand = c(0, 0)) +
  labs(x = sprintf("Genes colocalized with liver eQTL (SuSiE PP.H4 > 0.5); n = %d genes", n_tot), y = NULL,
       caption = sprintf("Finite SuSiE results in both classes (%d): %d enzyme-supported, %d direct-supported, %d both",
                         both_tested, both_enz, both_dir, enz_pos[["supported"]])) +
  guides(fill = guide_legend(ncol = 1), colour = guide_legend(ncol = 1)) +
  theme_masld(base_size = 6) +
  theme(axis.ticks.y = element_blank(), legend.position = "bottom", legend.justification = "left",
        legend.key.size = unit(0.22, "cm"), legend.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(size = 6, colour = "black", hjust = 0, face = "plain",
                                    margin = margin(t = 4)),
        plot.caption.position = "plot", legend.location = "plot")
ggsave(pdf_out, plot, width = 3.9, height = 2.1, device = cairo_pdf, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(out, "logs/sessionInfo_fig1c.txt"))
message(sprintf("Caption: %d genes colocalized with liver eQTL. ALT/AST/GGT (%d studies): %d genes = %d also direct-supported + %d direct-tested at PP.H4 <= 0.5 + %d not SuSiE-tested in a direct study; %d genes SuSiE-tested. NAFLD/NASH/PDFF (%d studies): %d = %d + %d + %d; %d tested. Tested in both: %d; %d enzyme-supported, %d direct-supported.",
  n_tot, n_enz_studies, sum(enz_pos), enz_pos[1], enz_pos[2], enz_pos[3], tested_enz, n_dir_studies,
  sum(dir_pos), dir_pos[1], dir_pos[2], dir_pos[3], tested_dir, both_tested, both_enz, both_dir))
cat("FIG1C_CANDIDATE\tPASS\n")
