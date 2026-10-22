#!/usr/bin/env Rscript
# F2abce candidate for Fig 2B (findings R02-04, R02-05). Derived from
# scripts/figures/fig2_resolution_plane.R (copy: fig2_resolution_plane.orig.R).
#
# Changes from the installed panel:
#  x  = within-ancestry fine-mapping PIP of the SAME variant that gives y:
#       the lead SNP (top_snp) of the gene's max-PP.H4.susie row, in that GWAS.
#       Within-ancestry PIP follows the aggregator's own rule before joint
#       integration (04_aggregate_results.R:172-184): CARMA for non-EUR strata
#       or unreliable SuSiE; max(SuSiE, CARMA) otherwise. SuSiEx / MESuSiE
#       cross-ancestry joint PIPs are not used. No gene-level max over GWAS.
#  y  = unchanged: max PP.H4.susie over the 35 main strata, COLOC r2.
#  colour = TxDb hg19 location class of that same lead variant, with coding
#       split into the plotted gene's own CDS/splice site versus another gene's.
#       No curated class overrides.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})
set.seed(20260923)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WS   <- file.path(BASE, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F2abce")
OUTF <- file.path(BASE, "figures/candidates/review-r3-fixes-20261015T112148Z/F2abce")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(WS, "scripts/lead_coding_owner.R"))
FM <- file.path(BASE, "GWAS/finemapping/results")
COLOC_INPUT <- file.path(FM, "susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv")

MAIN <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[placement == "main", study_name]
coloc <- fread(COLOC_INPUT, select = c("gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp"))
coloc <- coloc[gwas_name %in% MAIN & !is.na(gene) & gene != "" & !is.na(top_snp) & top_snp != ""]
setnames(coloc, c("PP.H4.susie", "PP.H4.abf"), c("pp4_susie", "pp4_abf"))
coloc[, pp4_susie := suppressWarnings(as.numeric(pp4_susie))]
coloc[, c("snp_chr", "snp_pos") := tstrsplit(top_snp, ":", fixed = TRUE)]
coloc[, snp_pos := as.integer(snp_pos)]
cat("main strata:", length(MAIN), " coloc rows:", nrow(coloc), " genes:", uniqueN(coloc$gene), "\n")

# y and its source row (same tie rule as the producer: first max in file order)
max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
at_max <- function(v, lab) if (all(is.na(v))) NA_character_ else lab[which.max(fifelse(is.na(v), -Inf, v))]
gy <- coloc[, .(pp4_susie = max_or_na(pp4_susie),
                y_gwas = at_max(pp4_susie, gwas_name),
                y_top_snp = at_max(pp4_susie, top_snp)), by = gene]
gy <- gy[!is.na(pp4_susie)]
gy[, c("ychr", "ypos") := tstrsplit(y_top_snp, ":", fixed = TRUE)]
gy[, ypos := as.integer(ypos)]
keys <- unique(gy[, .(study = y_gwas, chromosome = ychr, position = ypos)])
cat("genes with finite SuSiE PP.H4:", nrow(gy), " distinct y lead keys:", nrow(keys), "\n")

cols <- c("study", "chromosome", "position", "ancestry", "susie_pip", "carma_pip",
          "susie_reliable", "recommended_pip", "susiex_pip", "mesusie_pip")
canon <- fread(file.path(FM, "combined_finemapping.csv"), select = cols,
               colClasses = list(character = "chromosome"))
canon <- canon[keys, on = .(study, chromosome, position), nomatch = 0L]
canon[, x_source := "canonical_finemapping"]
targ <- fread(file.path(FM, "combined_finemapping_coloc_targeted.csv"),
              select = setdiff(cols, c("susiex_pip", "mesusie_pip")),
              colClasses = list(character = "chromosome"))
targ <- targ[keys, on = .(study, chromosome, position), nomatch = 0L]
targ[, `:=`(susiex_pip = NA_real_, mesusie_pip = NA_real_, x_source = "coloc_targeted_finemapping")]
fm <- rbindlist(list(canon, targ), use.names = TRUE)
num <- c("susie_pip", "carma_pip", "recommended_pip", "susiex_pip", "mesusie_pip")
fm[, (num) := lapply(.SD, function(v) suppressWarnings(as.numeric(v))), .SDcols = num]
fm[, susie_reliable := as.logical(susie_reliable)]
fm[, wa_pip := fifelse(ancestry != "EUR" | !susie_reliable %in% TRUE, carma_pip,
                       pmax(susie_pip, carma_pip, na.rm = TRUE))]
fm[, producer_pip := pmax(recommended_pip, susie_pip, na.rm = TRUE)]
fm[, source_rank := fifelse(x_source == "canonical_finemapping", 1L, 2L)]
# one row per key: largest within-ancestry PIP, canonical first on ties
fm <- fm[order(study, chromosome, position, -fifelse(is.na(wa_pip), -Inf, wa_pip), source_rank)]
fx <- fm[, .(x_wa = wa_pip[1], x_ancestry = ancestry[1], x_source = x_source[1],
             x_susie = max_or_na(susie_pip), x_carma = max_or_na(carma_pip),
             x_producer_rule = max_or_na(producer_pip)),
         by = .(study, chromosome, position)]
dt <- merge(gy, fx, by.x = c("y_gwas", "ychr", "ypos"), by.y = c("study", "chromosome", "position"),
            all.x = TRUE)

# Installed panel values for comparison
inst <- fread(file.path(BASE, "figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc_source.csv"))
dt <- merge(dt, inst[, .(gene, installed_x = causal_pip, installed_x_top_snp = x_top_snp,
                         installed_class = class, installed_plotted = plotted)],
            by = "gene", all.x = TRUE)

# Lead-variant class at the y lead
ann <- annotate_leads(dt$gene, dt$y_top_snp)
dt <- merge(dt, ann, by.x = c("gene", "y_top_snp"), by.y = c("gene", "variant_key"), all.x = TRUE)
dt[, class := fifelse(txdb_class %chin% c("coding", "spliceSite") & own_gene_coding,
                      "Coding, this gene",
              fifelse(txdb_class %chin% c("coding", "spliceSite"), "Coding, other gene",
                      "Non-coding"))]

roster <- c("HSD17B13", "TM6SF2", "MBOAT7", "GCKR", "RORA")
PP4_GATE <- 0.05; PIP_GATE <- 0.99
dt[, has_x := is.finite(x_wa)]
dt[, plotted := has_x & ((pp4_susie > PP4_GATE) | (x_wa >= PIP_GATE) | gene %in% roster)]
dt[, pp4 := pmax(pp4_susie, 0)]

# Summary numbers
s <- list(
  genes_finite_pp4 = nrow(dt),
  genes_with_x_same_row = dt[has_x == TRUE, .N],
  plotted = dt[plotted == TRUE, .N],
  installed_plotted = inst[plotted == TRUE, .N],
  corner_pip_ge_0.9_pp4_lt_0.5 = dt[plotted == TRUE & x_wa >= 0.9 & pp4 < 0.5, .N],
  installed_corner = inst[plotted == TRUE & causal_pip >= 0.9 & pp4 < 0.5, .N],
  topright_pip_ge_0.9_pp4_gt_0.5 = dt[plotted == TRUE & x_wa >= 0.9 & pp4 > 0.5, .N],
  installed_topright = inst[plotted == TRUE & causal_pip >= 0.9 & pp4 > 0.5, .N],
  plotted_noncoding = dt[plotted == TRUE & class == "Non-coding", .N],
  plotted_coding_this_gene = dt[plotted == TRUE & class == "Coding, this gene", .N],
  plotted_coding_other_gene = dt[plotted == TRUE & class == "Coding, other gene", .N],
  installed_class_counts = paste(inst[plotted == TRUE, .N, by = class][, paste0(class, "=", N)], collapse = ";"),
  susie_only_same_row_plotted = dt[is.finite(x_susie) & ((pp4 > PP4_GATE) | (x_susie >= PIP_GATE) | gene %in% roster), .N],
  susie_only_same_row_corner = dt[is.finite(x_susie) & ((pp4 > PP4_GATE) | (x_susie >= PIP_GATE) | gene %in% roster) & x_susie >= 0.9 & pp4 < 0.5, .N],
  producer_rule_same_row_plotted = dt[is.finite(x_producer_rule) & ((pp4 > PP4_GATE) | (x_producer_rule >= PIP_GATE) | gene %in% roster), .N],
  producer_rule_same_row_corner = dt[is.finite(x_producer_rule) & ((pp4 > PP4_GATE) | (x_producer_rule >= PIP_GATE) | gene %in% roster) & x_producer_rule >= 0.9 & pp4 < 0.5, .N]
)
for (n in names(s)) cat(sprintf("%-40s %s\n", n, s[[n]]))
cat("\nRoster:\n"); print(dt[gene %in% roster, .(gene, y_gwas, y_top_snp, pp4, x_wa, x_susie, x_carma, x_ancestry,
                                                  installed_x, installed_class, class, txdb_class, coding_symbols, plotted)])
cat("\nInstalled-plotted genes by x status under same-row rule:\n")
print(dt[installed_plotted == TRUE, .N, by = .(has_x, plotted)])

plot_dt <- dt[plotted == TRUE]
plot_dt[, class := factor(class, levels = c("Non-coding", "Coding, this gene", "Coding, other gene"))]
bg_dt <- plot_dt[!gene %in% roster][order(class)]
lab_dt <- plot_dt[gene %in% roster]
class_cols <- c("Non-coding" = "#1565C0", "Coding, this gene" = "#C9265E", "Coding, other gene" = "#7B1FA2")

p <- ggplot(plot_dt, aes(x = x_wa, y = pp4)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "grey45") +
  rasterize_layer(geom_point(data = bg_dt, aes(color = class), size = 0.45, alpha = 0.8, shape = 16)) +
  geom_point(data = lab_dt, aes(color = class), size = 1.6, shape = 16, alpha = 0.95, stroke = 0) +
  ggrepel::geom_text_repel(data = lab_dt, aes(label = gene), size = 6 / ggplot2::.pt,
    fontface = "plain", color = "black", max.overlaps = Inf, show.legend = FALSE,
    box.padding = 0.4, point.padding = 0.15, min.segment.length = 0, seed = 1,
    segment.size = 0.15, segment.color = "grey55", bg.color = "white", bg.r = 0.12) +
  scale_color_manual(values = class_cols, name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(-0.02, 1.03), breaks = c(0, 0.5, 1.0),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(-0.02, 1.03), breaks = c(0, 0.5, 1.0),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Within-ancestry PIP of COLOC lead", y = "Multi-signal COLOC PP.H4") +
  coord_cartesian(clip = "off") +
  theme_masld(base_size = 6) +
  theme(legend.position = "top", legend.direction = "vertical",
        legend.justification = "left",
        plot.margin = margin(2, 3, 1, 3, "pt"),
        axis.title = element_text(size = 6, face = "plain", color = "black"),
        axis.text = element_text(size = 6, face = "plain", color = "black"),
        legend.key.size = unit(5, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 6, face = "plain", color = "black")) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 0.95)))

pdf_file <- file.path(OUTF, "Fig2B_PIP_vs_SuSiE-coloc_sameRow.pdf")
ggsave(pdf_file, p, width = 1.80, height = 2.35, units = "in", device = cairo_pdf)

out <- dt[, .(gene, y_gwas, y_top_snp, pp4_susie = round(pp4_susie, 4),
              x_within_ancestry_pip = round(x_wa, 4), x_ancestry, x_source,
              x_susie_same_row = round(x_susie, 4), x_carma_same_row = round(x_carma, 4),
              x_installed_rule_same_row = round(x_producer_rule, 4),
              installed_x = installed_x, installed_x_top_snp, installed_class, installed_plotted,
              lead_txdb_class = txdb_class, lead_coding_genes = coding_symbols,
              own_gene_coding, class, is_label = gene %in% roster & plotted, plotted)]
setorder(out, -plotted, gene)
fwrite(out, file.path(OUTF, "Fig2B_PIP_vs_SuSiE-coloc_sameRow_source.csv"))
fwrite(data.table(metric = names(s), value = vapply(s, as.character, "")),
       file.path(OUTF, "Fig2B_sameRow_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(WS, "sessionInfo_fig2b.txt"))
cat("done\n")
