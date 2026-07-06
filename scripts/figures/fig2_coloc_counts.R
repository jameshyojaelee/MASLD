#!/usr/bin/env Rscript
# fig2_coloc_counts.R  (2026-06-12; cross-ancestry split 2026-07-05)  — Figure 2B
# GWAS-eQTL colocalization nominations: coloc.abf (baseline) vs SuSiE-coloc
# (primary) gene counts at PP.H4 thresholds, across the Broadaway cis-eQTL genes.
# Scoped 2026-07-06 to the 35 Tier-1/2 (liver-specific) MAIN strata only (NAFLD/NASH/
# PDFF + ALT/AST/GGT; placement=="main" in gwas_trait_tier.tsv); Tier-3/4 supp strata
# move to a supplementary full-portfolio figure. The methods cross over (ABF
# nominates more at > 0.5, SuSiE-coloc is more confident at > 0.9).
# Each method's count is split into EUR-headline vs cross-ancestry-headline
# (headline PP.H4 driven by a NON-EUR GWAS; method-matched flags coloc_abf/
# coloc_susie_headline_cross_anc). The eQTL panel is EUR, so a non-EUR-GWAS
# headline holds a lower evidentiary bar — drawn as a faded cap on each bar.
# Counts computed from disk so they always match the manuscript text (never hardcoded).
#
# Out: figures/main/fig2_genetics/panels/coloc_method_counts.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

# ── MAIN (Tier-1/2, liver-specific) restriction (2026-07-06) ──────────────────
# Scoped to the Tier-1/2 MAIN strata (placement=="main" in the tier map); the Tier-3/4
# supp strata (MVP Cirrhosis/ChronLiver/Albumin/Platelet) move to a supplementary
# full-portfolio figure. The pre-aggregated wide gene_level_coloc.csv carries
# FULL-portfolio best columns, so we recompute each gene's best-SuSiE / best-ABF PP.H4
# (and the ancestry that drives that best) over the MAIN strata directly from the long
# per-gene-per-GWAS coloc table.
MAIN_STUDIES <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[
  placement == "main", study_name]
# cis-eQTL gene universe (denominator) is a property of the Broadaway EUR eQTL panel,
# independent of the GWAS tiering — read its row count once.
n_genes <- nrow(fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"),
                      select = "gene"))
sc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
            select = c("gwas_name", "gene", "PP.H4.susie", "PP.H4.abf"))
sc <- sc[gwas_name %in% MAIN_STUDIES]
sc[, ancestry := as.character(gwas_ancestry(gwas_name))]
sc[, `:=`(PP.H4.abf   = suppressWarnings(as.numeric(PP.H4.abf)),
          PP.H4.susie = suppressWarnings(as.numeric(PP.H4.susie)))]
# per-gene best ABF / best SuSiE over MAIN strata + the ancestry that drives that best
abf_g <- sc[!is.na(PP.H4.abf)][order(gene, -PP.H4.abf), .SD[1], by = gene][
  , .(gene, abf = PP.H4.abf, abf_anc = ancestry)]
su_g  <- sc[!is.na(PP.H4.susie)][order(gene, -PP.H4.susie), .SD[1], by = gene][
  , .(gene, su = PP.H4.susie, su_anc = ancestry)]
d   <- merge(abf_g, su_g, by = "gene", all = TRUE)
abf <- d$abf; su <- d$su
# Cross-ancestry provenance: headline PP.H4 driven by a NON-EUR GWAS (method-matched).
# The eQTL panel is EUR, so a non-EUR-GWAS headline holds a lower evidentiary bar — split
# each method's colocalizing-gene count into EUR-headline vs cross-anc.
abf_xa <- !is.na(d$abf_anc) & d$abf_anc != "EUR"
su_xa  <- !is.na(d$su_anc)  & d$su_anc  != "EUR"

ths <- c(0.5, 0.8, 0.9)
mk <- function(t) data.table(
  threshold  = sprintf("PP.H4 > %.1f", t),
  method     = c("coloc.abf (baseline)", "coloc.abf (baseline)",
                 "SuSiE-coloc (primary)", "SuSiE-coloc (primary)"),
  provenance = c("EUR", "cross-ancestry", "EUR", "cross-ancestry"),
  n          = c(sum(abf > t & !abf_xa, na.rm = TRUE), sum(abf > t & abf_xa, na.rm = TRUE),
                 sum(su  > t & !su_xa,  na.rm = TRUE), sum(su  > t & su_xa,  na.rm = TRUE)))
counts <- rbindlist(lapply(ths, mk))
counts[, threshold := factor(threshold, levels = sprintf("PP.H4 > %.1f", ths))]
counts[, method := factor(method, levels = c("coloc.abf (baseline)", "SuSiE-coloc (primary)"))]
# EUR at the base, cross-ancestry stacked on top (first level draws on top)
counts[, provenance := factor(provenance, levels = c("cross-ancestry", "EUR"))]
# grouped (dodge by method) + stacked (by provenance): numeric x offset per method
counts[, xpos := as.integer(threshold) + ifelse(method == "coloc.abf (baseline)", -0.19, 0.19)]
totals <- counts[, .(n = sum(n)), by = .(threshold, method, xpos)]  # bar-top labels

method_cols <- c("coloc.abf (baseline)" = "#90A4AE", "SuSiE-coloc (primary)" = "#1565C0")

p <- ggplot(counts, aes(x = xpos, y = n, fill = method, alpha = provenance)) +
  geom_col(width = 0.34) +
  geom_text(data = totals, aes(x = xpos, y = n, label = n), inherit.aes = FALSE,
            vjust = -0.35, size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = method_cols, name = NULL) +
  scale_alpha_manual(values = c("EUR" = 1, "cross-ancestry" = 0.42),
                     breaks = c("EUR", "cross-ancestry"),
                     labels = c("EUR" = "EUR", "cross-ancestry" = "cross-ancestry (non-EUR GWAS)"),
                     name = "COLOC headline") +
  scale_x_continuous(breaks = seq_along(levels(counts$threshold)),
                     labels = levels(counts$threshold)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Colocalizing genes",
       title = "GWAS-eQTL colocalization (Broadaway liver eQTLs)",
       subtitle = sprintf("%s cis-eQTL genes, 35 Tier-1/2 (liver-specific) GWAS; bar height = colocalizing genes, faded cap = cross-ancestry headline",
                          format(n_genes, big.mark = ","))) +
  theme_masld(base_size = 9) +
  theme(legend.position = c(0.98, 0.97), legend.justification = c(1, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
        plot.title = element_text(size = 9),
        plot.subtitle = element_text(size = 6.5, color = "grey30"))

save_fig(p, file.path(PANEL_DIR, "coloc_method_counts.pdf"),
         width = fig_col_width * 1.05, height = 2.9)

fwrite(dcast(counts, threshold ~ method + provenance, value.var = "n"),
       file.path(PANEL_DIR, "coloc_method_counts_source.csv"))
message("CAPTION: GWAS-eQTL colocalization across the 35 Tier-1/2 (liver-specific) GWAS ",
        "(Broadaway EUR liver cis-eQTLs). coloc.abf (baseline) vs SuSiE-coloc (primary); ",
        "bar height = colocalizing genes at each PP.H4 threshold, split into EUR-headline ",
        "(solid) vs cross-ancestry-headline (faded, headline PP.H4 from a non-EUR GWAS). ",
        "Because the eQTL panel is EUR, cross-ancestry headlines hold a lower evidentiary bar.")
cat(sprintf("[fig2B] wrote coloc_method_counts.pdf (n=%d genes; PP.H4>0.5: ABF %d [%d cross-anc] vs SuSiE %d [%d cross-anc])\n",
            n_genes,
            sum(abf > 0.5, na.rm = TRUE), sum(abf > 0.5 & abf_xa, na.rm = TRUE),
            sum(su  > 0.5, na.rm = TRUE), sum(su  > 0.5 & su_xa,  na.rm = TRUE)))
