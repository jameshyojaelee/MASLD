#!/usr/bin/env Rscript
# F2abce repair render of the Fig 2B candidate (R02-04, R02-05).
# Same data rule as the round-1 candidate (x = largest SuSiE PIP at y's GWAS and
# lead position; y = max PP.H4.susie over the 35 main strata, COLOC r2).
# Changes from round 1: x re-derived from the release files in job 7520
# (fig2b_multirow_check.R); lead class re-annotated here; the class key prints
# the gene count of each class and the total; points drawn as vectors.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(20260923)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RD   <- file.path(BASE, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce")
OUTF <- file.path(BASE, "figures/candidates/review-r3-fixes-20261015T112148Z/F2abce/repair")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(RD, "scripts/lead_coding_owner.R"))

dt <- fread(file.path(RD, "out/fig2b_same_row_x_by_rule.tsv"), colClasses = list(character = "chromosome"))
roster <- c("HSD17B13", "TM6SF2", "MBOAT7", "GCKR", "RORA")
pl <- dt[plotted_max == TRUE]
ann <- annotate_leads(pl$gene, pl$y_top_snp)
pl <- merge(pl, ann, by.x = c("gene", "y_top_snp"), by.y = c("gene", "variant_key"), all.x = TRUE)
pl[, class := fifelse(txdb_class %chin% c("coding", "spliceSite") & own_gene_coding, "coding, this gene",
              fifelse(txdb_class %chin% c("coding", "spliceSite"), "coding, other gene", "noncoding"))]
pl[, pp4p := pmax(pp4, 0)]
lv <- c("noncoding", "coding, this gene", "coding, other gene")
n_by <- pl[, .N, by = class]
labs_cls <- setNames(paste0(lv, " (", n_by[match(lv, class), N], ")"), lv)
pl[, class := factor(class, levels = lv)]
cat("plotted:", nrow(pl), "\n"); print(n_by)
cat("corner PIP>=0.9 & PP.H4<0.5:", pl[x_max >= 0.9 & pp4p < 0.5, .N],
    " top-right PIP>=0.9 & PP.H4>0.5:", pl[x_max >= 0.9 & pp4p > 0.5, .N], "\n")
print(pl[class != "noncoding", .(gene, class, y_top_snp, txdb_class, coding_symbols)][order(class, gene)])

bg_dt <- pl[!gene %in% roster][order(class)]
lab_dt <- pl[gene %in% roster]
class_cols <- c("noncoding" = "#1565C0", "coding, this gene" = "#C9265E", "coding, other gene" = "#7B1FA2")
p <- ggplot(pl, aes(x = x_max, y = pp4p)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "grey45") +
  geom_point(data = bg_dt, aes(color = class), size = 0.45, alpha = 0.8, shape = 16, stroke = 0) +
  geom_point(data = lab_dt, aes(color = class), size = 1.6, shape = 16, alpha = 0.95, stroke = 0) +
  ggrepel::geom_text_repel(data = lab_dt, aes(label = gene), size = 6 / ggplot2::.pt,
    fontface = "plain", color = "black", max.overlaps = Inf, show.legend = FALSE,
    box.padding = 0.4, point.padding = 0.15, min.segment.length = 0, seed = 1,
    segment.size = 0.15, segment.color = "grey55") +
  scale_color_manual(values = class_cols, labels = labs_cls, name = paste0(nrow(pl), " genes"), drop = FALSE) +
  scale_x_continuous(limits = c(-0.02, 1.03), breaks = c(0, 0.5, 1.0), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(-0.02, 1.03), breaks = c(0, 0.5, 1.0), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "SuSiE PIP of COLOC lead", y = "Multi-signal COLOC PP.H4") +
  coord_cartesian(clip = "off") +
  theme_masld(base_size = 6) +
  theme(legend.position = "top", legend.direction = "vertical", legend.justification = "left",
        plot.margin = margin(2, 3, 1, 3, "pt"),
        axis.title = element_text(size = 6, face = "plain", color = "black"),
        axis.text = element_text(size = 6, face = "plain", color = "black"),
        legend.title = element_text(size = 6, face = "plain", color = "black"),
        legend.key.size = unit(5, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 6, face = "plain", color = "black")) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 0.95)))
pdf_file <- file.path(OUTF, "Fig2B_PIP_vs_SuSiE-coloc_sameRow_susie_repair.pdf")
ggsave(pdf_file, p, width = 1.80, height = 2.45, units = "in", device = cairo_pdf)
out <- pl[, .(gene, y_gwas, y_top_snp, pp4_susie = round(pp4, 4), x_susie_max = round(x_max, 4),
              x_susie_min = round(x_min, 4), n_finemapping_rows = n_rows, n_distinct_pip,
              lead_txdb_class = txdb_class, lead_coding_genes = coding_symbols, own_gene_coding,
              class, is_label = gene %in% roster)]
setorder(out, gene)
fwrite(out, file.path(OUTF, "Fig2B_PIP_vs_SuSiE-coloc_sameRow_susie_repair_source.csv"))
writeLines(capture.output(sessionInfo()), file.path(RD, "out/sessionInfo_fig2b_render.txt"))
cat("done\n")
