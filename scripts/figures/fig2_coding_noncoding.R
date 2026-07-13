#!/usr/bin/env Rscript
# fig2_coding_noncoding.R  (2026-06-12)  — Fig 2C
# Where does the colocalizing causal variant fall? Coding vs non-coding split.
#
# Answers the draft's Fig 2C placeholder "[How many noncoding vs coding variants?]".
# Single horizontal stacked bar over the de-duplicated SuSiE-OR-ABF colocalization set
# (best PP.H4 > 0.5 per gene), segmented by the VEP CONSEQUENCE of that gene's strongest
# COLOC lead SNP (NOT a TSS-distance cut). Headline (MAIN, Tier-1/2): only ~34/1,031 (~3%)
# colocalizing genes are coding-led; ~97% act through non-coding regulatory sequence ->
# motivates the expression-mediated effector model (coding class deferred to later figs).
#
# 2026-06-25: switched from SuSiE-only to the de-duped SuSiE-OR-ABF union so the
#             denominator MATCHES Fig2G (ancestry specificity, also SuSiE/ABF). Per-gene
#             consequence = the fine_class of the gene's max-PP.H4 colocalization.
# 2026-07-04: portfolio rebuilt with MVP (23->50 GWAS).
# 2026-07-06: scoped to the Tier-1/2 MAIN strata only (placement=="main"); union drops from
#             the full-portfolio 1,527 to the MAIN 1,031 (coding 59->34). All counts recomputed
#             from disk, so the plot tracks the restriction automatically.
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
# MAIN (Tier-1/2, liver-specific) restriction (2026-07-06): keep only variant/gene rows
# whose driving `study` is a placement=="main" stratum (NAFLD/NASH/PDFF + ALT/AST/GGT);
# Tier-3/4 supp strata move to a supplementary full-portfolio figure. Keeps the denominator
# aligned with the MAIN-restricted Fig2G ancestry-specificity panel.
MAIN_STUDIES <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[
  placement == "main", study_name]
av <- av[study %in% MAIN_STUDIES]
setorder(av, gene_symbol, -pp4_best)
best <- av[, .SD[1], by = gene_symbol][pp4_best > 0.5]
# Merge the single spliceSite gene (1 gene) into coding to simplify the legend and bar
best[fine_class == "spliceSite", fine_class := "coding"]
best[fine_class %in% c("fiveUTR", "threeUTR"), fine_class := "UTR"]   # clump 5'/3' UTR into one group
f1 <- best[, .(n_genes = .N), by = fine_class]
f1[, `:=`(method = "SuSiE or ABF",                 # single bar over the union set
          coarse_class = fifelse(fine_class == "coding", "coding", "non-coding"))]

# consequence order — coding first so that coding lands at the very top of the vertical stack
lev <- c("coding","UTR","promoter","intron","intergenic")
labs <- c(intergenic="intergenic", intron="intron", promoter="promoter",
          UTR="UTR", coding="coding")
f1[, fine_class := factor(fine_class, levels = lev)]
f1[, method := factor(method, levels = "SuSiE or ABF")]
present_lev <- lev[lev %in% as.character(f1$fine_class)]   # legend shows only classes with data

cls_cols <- c(intergenic="#B0BEC5", intron="#64B5F6", promoter="#1565C0",
              UTR="#4DB6AC", coding="#C9265E")

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
            color = "white", fontface = "plain", size = GEOM_TEXT_6PT) +
  # total above the vertical bar (black)
  geom_text(data = tot, aes(x = method, y = n, label = n),
            inherit.aes = FALSE, vjust = -0.6, hjust = 0.5, size = GEOM_TEXT_6PT, fontface = "plain") +
  # coding/non-coding bracket labels dropped for the slim 0.97in-wide column (no room
  # for side labels); the split (non-coding 997 / coding 34) is in the caption. Visually:
  # coding = the magenta sliver at top, everything below = non-coding.
  scale_fill_manual(values = cls_cols, labels = unname(labs[present_lev]),
                    name = NULL, breaks = present_lev) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  scale_x_discrete(expand = expansion(add = c(0.35, 0.35))) +
  labs(x = NULL, y = "Colocalizing genes") +   # PP.H4>0.5 threshold in caption (title too long for the short vertical axis)
  coord_cartesian(clip = "off") +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom", legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = 6), legend.margin = margin(1, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"), legend.spacing.x = unit(1, "pt"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 1))   # 0.96in width fits only a single-column legend (2-col words overflow)

message(sprintf("[caption] Fig2C variant class: non-coding %d (%.0f%%) vs coding %d (%.0f%%); best PP.H4 > 0.5.",
                nnon, 100 - ann$pct, ncod, ann$pct))
source(file.path(BASE, "figures/layout_specs/regenerate_panels.R"))   # save_panel(): exact contract size + cairo_pdf
save_panel(p, "main/fig2_genetics/panels/Fig2C_coding_noncoding_split.pdf",
           read_sizes(file.path(BASE, "figures/layout_specs/figure2_panel_sizes.tsv")),
           file.path(BASE, "figures"))

out <- merge(f1[, .(method, consequence = as.character(fine_class), coarse_class, n_genes)],
             ann[, .(method, method_total = n, coding_led = coding, coding_pct = pct)],
             by = "method")
fwrite(out, file.path(PANEL_DIR, "Fig2C_coding_noncoding_source.csv"))
cat(sprintf("[fig2C] wrote Fig2C_coding_noncoding_split.pdf  (SuSiE-or-ABF coding-led = %d/%d = %.1f%%)\n",
            ann$coding, ann$n, ann$pct))
message(sprintf("CAPTION (Fig2C): Consequence class of each colocalizing gene's strongest COLOC lead SNP (best PP.H4 > 0.5, SuSiE or ABF; de-duplicated; same %d-gene set as Fig2G). Only %d (%.1f%%) are coding-led; %.0f%% act through non-coding regulatory sequence.",
                ann$n, ann$coding, ann$pct, 100 - ann$pct))
