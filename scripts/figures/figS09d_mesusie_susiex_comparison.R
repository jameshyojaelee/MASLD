#!/usr/bin/env Rscript
# figS09d_mesusie_susiex_comparison.R
# Compact 4-panel supplementary figure comparing MESuSiE and SuSiEX
# cross-ancestry fine-mapping results (EUR UKBB × EAS BBJ, ALT/AST/GGT)
#
# Panel A: Per-trait locus architecture — shared/EUR-only/EAS-only CS (MESuSiE)
# Panel B: Variant-level PIP scatter (MESuSiE vs SuSiEX)
# Panel C: Credible set size comparison — MESuSiE shared CS vs SuSiEX joint CS
# Panel D: Known MASLD gene concordance heatmap
#
# Usage: Rscript scripts/figures/figS09d_mesusie_susiex_comparison.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FM_DIR   <- file.path(BASE, "GWAS/finemapping/results")
OUT_DIR  <- FIGS09_DIR
dir.create(file.path(OUT_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# Ancestry colors
eur_col    <- "#1565C0"   # Blue
eas_col    <- "#C2185B"   # Magenta
shared_col <- "#00695C"   # Teal (shared EUR+EAS)
mesusie_col <- "#7B1FA2"  # Violet
susiex_col  <- "#E65100"  # Deep orange

cat("Loading data...\n")

# ---------------------------------------------------------------------------
# Load MESuSiE data
# ---------------------------------------------------------------------------
me_locus <- fread(file.path(FM_DIR, "mesusie/mesusie_locus_summary.csv"))
me_gene  <- fread(file.path(FM_DIR, "mesusie/mesusie_gene_summary.csv"))
me_var   <- fread(file.path(FM_DIR, "mesusie/mesusie_variant_summary.csv"))

# ---------------------------------------------------------------------------
# Load SuSiEX data
# ---------------------------------------------------------------------------
sx_loci  <- fread(file.path(FM_DIR, "susiex/shared_loci.csv"))
sx_gene  <- fread(file.path(FM_DIR, "susiex/susiex_gene_summary.csv"))
sx_var   <- fread(file.path(FM_DIR, "susiex/susiex_variant_summary.csv"))

cat(sprintf("MESuSiE: %d loci, %d genes, %d variants\n",
            nrow(me_locus), nrow(me_gene), nrow(me_var)))
cat(sprintf("SuSiEX:  %d loci, %d genes, %d variants\n",
            nrow(sx_loci), nrow(sx_gene), nrow(sx_var)))

# ===========================================================================
# Panel A: Credible set architecture per trait (MESuSiE)
# Stacked bar: shared / EUR-only / EAS-only CS counts per trait
# ===========================================================================
cat("Panel A: CS architecture...\n")

cs_summary <- me_locus[, .(
  Shared   = sum(n_cs_shared, na.rm = TRUE),
  EUR_only = sum(n_cs_eur, na.rm = TRUE),
  EAS_only = sum(n_cs_eas, na.rm = TRUE),
  n_loci   = .N
), by = trait_pair]

cs_long <- melt(cs_summary, id.vars = c("trait_pair", "n_loci"),
                variable.name = "CS_type", value.name = "n_cs")
cs_long[, CS_type := factor(CS_type, levels = c("EAS_only", "EUR_only", "Shared"))]
cs_long[, trait_label := paste0(trait_pair, "\n(", n_loci, " loci)")]

pA <- ggplot(cs_long, aes(x = trait_label, y = n_cs, fill = CS_type)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(Shared = shared_col, EUR_only = eur_col, EAS_only = eas_col),
    labels = c(Shared = "Shared (EUR+EAS)", EUR_only = "EUR-specific",
               EAS_only = "EAS-specific")
  ) +
  labs(x = NULL, y = "Number of credible sets",
       fill = NULL) +
  theme_masld() +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6)
  )
message("[caption] Panel A: Credible set architecture (MESuSiE)")

# ===========================================================================
# Panel B: Variant PIP scatter (MESuSiE vs SuSiEX)
# ===========================================================================
cat("Panel B: PIP scatter...\n")

# Merge on locus_id + variant_id
me_var[, key := paste0(locus_id, "_", variant_id)]
sx_var[, key := paste0(locus_id, "_", variant_id)]

shared_vars <- merge(
  me_var[, .(key, mesusie_pip = pip, cs_cat = cs_category, in_cs_me = in_cs)],
  sx_var[, .(key, susiex_pip = OVRL_PIP)],
  by = "key"
)

cat(sprintf("Shared variant-locus pairs: %d\n", nrow(shared_vars)))

# Compute correlation stats
r_all <- cor(shared_vars$mesusie_pip, shared_vars$susiex_pip,
             method = "spearman", use = "complete.obs")
high <- shared_vars[mesusie_pip > 0.1 & susiex_pip > 0.1]
r_high <- if(nrow(high) > 5) cor(high$mesusie_pip, high$susiex_pip,
                                   method = "pearson") else NA

cat(sprintf("Spearman rho (all): %.3f\n", r_all))
cat(sprintf("Pearson r (both PIP>0.1, n=%d): %.3f\n", nrow(high), r_high))

# Color by MESuSiE CS category
shared_vars[, cs_label := fcase(
  cs_cat == "EUR_EAS", "Shared CS",
  cs_cat == "EUR",     "EUR-only CS",
  cs_cat == "EAS",     "EAS-only CS",
  default = "Not in CS"
)]
shared_vars[, cs_label := factor(cs_label,
  levels = c("Not in CS", "EUR-only CS", "EAS-only CS", "Shared CS"))]

# Subsample background points for plotting speed
set.seed(42)
bg <- shared_vars[cs_label == "Not in CS"]
bg_sub <- bg[sample(.N, min(.N, 5000))]
fg <- shared_vars[cs_label != "Not in CS"]
plot_dat <- rbind(bg_sub, fg)

pB <- ggplot(plot_dat, aes(x = mesusie_pip, y = susiex_pip, color = cs_label)) +
  geom_point(data = plot_dat[cs_label == "Not in CS"], size = 0.3, alpha = 0.15) +
  geom_point(data = plot_dat[cs_label != "Not in CS"], size = 0.8, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c(
    "Not in CS"    = "gray70",
    "EUR-only CS"  = eur_col,
    "EAS-only CS"  = eas_col,
    "Shared CS"    = shared_col
  )) +
  annotate("text", x = 0.05, y = 0.92,
           label = sprintf("rho==%.2f~~(all)", r_all),
           parse = TRUE, size = GEOM_TEXT_6PT, hjust = 0, color = "black") +
  annotate("text", x = 0.05, y = 0.82,
           label = sprintf("italic(r)==%.2f~~(PIP>0.1)", r_high),
           parse = TRUE, size = GEOM_TEXT_6PT, hjust = 0, color = shared_col) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "MESuSiE PIP", y = "SuSiEX PIP",
       color = NULL) +
  theme_masld() +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
message("[caption] Panel B: Variant PIP concordance")

# ===========================================================================
# Panel C: CS size distribution — MESuSiE shared CS vs trait
# ===========================================================================
cat("Panel C: CS size distribution...\n")

# CS sizes from locus summary (min_cs_size_shared is the smallest shared CS at each locus)
cs_sizes <- me_locus[n_cs_shared > 0, .(
  trait_pair,
  shared_cs_size = min_cs_size_shared,
  locus_id
)]
cs_sizes <- cs_sizes[!is.na(shared_cs_size)]

# Also get total number of CS per category
total_cs <- data.table(
  category = c("Shared\n(EUR+EAS)", "EUR-\nspecific", "EAS-\nspecific"),
  n_cs = c(sum(me_locus$n_cs_shared, na.rm = TRUE),
           sum(me_locus$n_cs_eur, na.rm = TRUE),
           sum(me_locus$n_cs_eas, na.rm = TRUE)),
  n_variants = c(nrow(me_var[cs_category == "EUR_EAS" & in_cs == TRUE]),
                 nrow(me_var[cs_category == "EUR" & in_cs == TRUE]),
                 nrow(me_var[cs_category == "EAS" & in_cs == TRUE]))
)
total_cs[, category := factor(category,
  levels = c("Shared\n(EUR+EAS)", "EUR-\nspecific", "EAS-\nspecific"))]
total_cs[, avg_size := n_variants / n_cs]

pC <- ggplot(total_cs, aes(x = category, y = n_cs,
                            fill = category)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d CS\n(%.0f vars)", n_cs, n_variants)),
            vjust = -0.3, size = GEOM_TEXT_6PT) +
  scale_fill_manual(values = c(
    "Shared\n(EUR+EAS)" = shared_col,
    "EUR-\nspecific" = eur_col,
    "EAS-\nspecific" = eas_col
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "Number of credible sets") +
  theme_masld() +
  theme(
    legend.position = "none"
  )
message("[caption] Panel C: Cross-ancestry signal decomposition")

# ===========================================================================
# Panel D: Known MASLD gene concordance
# ===========================================================================
cat("Panel D: Gene concordance...\n")

# Known MASLD genes to check
known_genes <- c("PNPLA3", "GCKR", "ALDH2", "HSD17B13", "TM6SF2",
                 "MBOAT7", "MARC1", "SERPINA1", "FADS1", "FADS2",
                 "ABCB4", "SAMM50", "ERLIN1", "ATXN2")

gene_compare <- data.table(gene = known_genes)

# Merge MESuSiE info
me_sub <- me_gene[GeneSymbol %in% known_genes,
                  .(gene = GeneSymbol,
                    mesusie_pip = mesusie_max_pip,
                    mesusie_shared = mesusie_in_shared_cs,
                    mesusie_traits = trait_pairs)]
gene_compare <- merge(gene_compare, me_sub, by = "gene", all.x = TRUE)

# Merge SuSiEX info
sx_sub <- sx_gene[GeneSymbol %in% known_genes,
                  .(gene = GeneSymbol,
                    susiex_pip = susiex_max_pip,
                    susiex_traits = trait_pairs)]
gene_compare <- merge(gene_compare, sx_sub, by = "gene", all.x = TRUE)

# Replace NAs
gene_compare[is.na(mesusie_pip), mesusie_pip := 0]
gene_compare[is.na(susiex_pip), susiex_pip := 0]
gene_compare[is.na(mesusie_shared), mesusie_shared := FALSE]

# Create long format for heatmap
gene_long <- rbind(
  gene_compare[, .(gene, method = "MESuSiE", pip = mesusie_pip)],
  gene_compare[, .(gene, method = "SuSiEX", pip = susiex_pip)]
)

# Order genes by max PIP across methods
gene_order <- gene_compare[, .(max_pip = pmax(mesusie_pip, susiex_pip, na.rm = TRUE)),
                            by = gene][order(-max_pip)]$gene
gene_long[, gene := factor(gene, levels = rev(gene_order))]
gene_long[, method := factor(method, levels = c("MESuSiE", "SuSiEX"))]

# Shared CS indicator
shared_indicator <- gene_compare[mesusie_shared == TRUE, .(gene)]
gene_long[, is_shared := gene %in% shared_indicator$gene]

pD <- ggplot(gene_long, aes(x = method, y = gene, fill = pip)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(pip > 0, sprintf("%.2f", pip), "")),
            size = GEOM_TEXT_6PT, color = ifelse(gene_long$pip > 0.6, "white", "black")) +
  # Mark shared CS genes
  geom_point(data = gene_long[is_shared == TRUE & method == "MESuSiE"],
             aes(x = method, y = gene),
             shape = 8, size = 1.5, color = shared_col, inherit.aes = FALSE) +
  scale_fill_gradient2(
    low = "white", mid = "#F48FB1", high = "#880E4F",
    midpoint = 0.5, limits = c(0, 1),
    name = "Max PIP"
  ) +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    axis.text.y = element_text(face = "italic", size = 6),
    panel.grid = element_blank()
  )
message("[caption] Panel D: Known MASLD gene fine-mapping")

# ===========================================================================
# Assemble composite figure
# ===========================================================================
cat("Assembling composite figure...\n")

composite <- (pA | pB) / (pC | pD) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "plain", size = 11))
  )

# Save composite
pdf_path <- file.path(OUT_DIR, "panels", "figS09d_mesusie_susiex_comparison.pdf")
ggsave(pdf_path, composite, width = fig_full_width, height = 8 * fig_full_width / 10, device = cairo_pdf)
cat(sprintf("Saved: %s\n", pdf_path))

# ===========================================================================
# Print summary statistics
# ===========================================================================
cat("\n════════════════════════════════════════════════════════\n")
cat("CROSS-ANCESTRY FINE-MAPPING COMPARISON SUMMARY\n")
cat("════════════════════════════════════════════════════════\n\n")

cat(sprintf("Shared loci analyzed:\n"))
cat(sprintf("  SuSiEX: %d loci (ALT: %d, AST: %d, GGT: %d)\n",
            nrow(sx_loci),
            nrow(sx_loci[trait_pair == "ALT"]),
            nrow(sx_loci[trait_pair == "AST"]),
            nrow(sx_loci[trait_pair == "GGT"])))
cat(sprintf("  MESuSiE: %d loci (11 dropped due to <50 aligned SNPs)\n\n",
            nrow(me_locus)))

cat(sprintf("MESuSiE credible sets:\n"))
cat(sprintf("  Shared (EUR+EAS): %d CS in %d loci (%.0f%% loci have shared signal)\n",
            sum(me_locus$n_cs_shared, na.rm = TRUE),
            sum(me_locus$n_cs_shared > 0, na.rm = TRUE),
            100 * mean(me_locus$n_cs_shared > 0, na.rm = TRUE)))
cat(sprintf("  EUR-specific: %d CS in %d loci\n",
            sum(me_locus$n_cs_eur, na.rm = TRUE),
            sum(me_locus$n_cs_eur > 0, na.rm = TRUE)))
cat(sprintf("  EAS-specific: %d CS in %d loci\n\n",
            sum(me_locus$n_cs_eas, na.rm = TRUE),
            sum(me_locus$n_cs_eas > 0, na.rm = TRUE)))

cat(sprintf("Variant PIP concordance (n=%d shared variants):\n", nrow(shared_vars)))
cat(sprintf("  Spearman rho (all variants): %.3f\n", r_all))
cat(sprintf("  Pearson r (both PIP>0.1, n=%d): %.3f\n", nrow(high), r_high))
n_agree_high <- nrow(shared_vars[mesusie_pip > 0.5 & susiex_pip > 0.5])
n_me_high    <- nrow(shared_vars[mesusie_pip > 0.5])
n_sx_high    <- nrow(shared_vars[susiex_pip > 0.5])
cat(sprintf("  Both PIP>0.5: %d | MESuSiE-only: %d | SuSiEX-only: %d\n\n",
            n_agree_high, n_me_high - n_agree_high, n_sx_high - n_agree_high))

cat(sprintf("Gene-level overlap:\n"))
me_genes <- unique(me_gene$GeneSymbol)
sx_genes <- unique(sx_gene$GeneSymbol)
shared_genes <- intersect(me_genes, sx_genes)
cat(sprintf("  MESuSiE: %d genes | SuSiEX: %d genes | Shared: %d (%.0f%%)\n",
            length(me_genes), length(sx_genes), length(shared_genes),
            100 * length(shared_genes) / length(union(me_genes, sx_genes))))

cat("\nKnown MASLD gene fine-mapping:\n")
print(gene_compare[, .(gene, mesusie_pip = round(mesusie_pip, 3),
                        susiex_pip = round(susiex_pip, 3),
                        shared_CS = mesusie_shared)])

cat("\nDone.\n")
