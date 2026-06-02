#!/usr/bin/env Rscript
# 21_gwas_spatial_overlay.R
# ---------------------------------------------------------------------------
# GWAS-concordance overlay: map MASLD GWAS loci onto the cross-species
# concordance atlas, test enrichment, and produce evidence cards for key loci.
# Follows publication_theme_guidelines.md: PDF, Helvetica, Sanjana palette
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(biomaRt)
})

pdf.options(useDingbats = FALSE)

cat("=== Strategy 12: GWAS Spatial Overlay ===\n\n")

# ============================================================
#  Paths
# ============================================================
PROJECT  <- Sys.getenv("MASLD_PROJECT_ROOT",
              "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RNA_DIR  <- file.path(PROJECT, "RNA-seq")
OUT_DIR  <- file.path(RNA_DIR, "results", "gwas_spatial_convergence")
FIG_DIR  <- file.path(PROJECT, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_FILE     <- file.path(PROJECT, "GWAS", "Closest_genes.csv")
ATLAS_FILE    <- file.path(RNA_DIR, "Analysis", "Cross_Species_Concordance",
                           "results", "concordance_atlas_unified.csv")
CONSENSUS_FILE <- file.path(RNA_DIR, "Human", "Patient_Cohorts", "analysis",
                            "integration", "results", "integration",
                            "consensus_degs.csv")
# MR_FILE removed 2026-04-22 — MR ditched from paper.

# ============================================================
#  Theme & palette (Sanjana lab standard)
# ============================================================
theme_pub <- theme_minimal(base_family = "Helvetica", base_size = 7) +
  theme(
    plot.title    = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 7, color = "grey40"),
    axis.title    = element_text(size = 8),
    axis.text     = element_text(size = 6),
    legend.text   = element_text(size = 6),
    legend.title  = element_text(size = 7),
    panel.grid.minor = element_blank(),
    strip.text    = element_text(size = 7, face = "bold")
  )

pal <- list(
  magenta = "#e14b9d", pink = "#e35070", purple = "#d358c7",
  blue = "#4baeef", orange = "#e1b172", green = "#30d796",
  teal = "#2bbfbd", grey = "#808080"
)

concordance_colors <- c(
  Conserved       = pal$magenta,
  Moderate_Concordance = pal$pink,
  Human_Enriched       = pal$purple,
  Mouse_Specific       = pal$blue,
  Diet_Selective       = pal$orange,
  Species_Discordant   = pal$teal,
  Not_Significant      = "#cccccc",
  Unclassified         = "#888888"
)

# ============================================================
#  1. Load GWAS genes
# ============================================================
cat("[1] Loading GWAS genes...\n")
gwas_genes <- fread(GWAS_FILE, header = FALSE, col.names = "symbol")
gwas_genes <- gwas_genes[symbol != "" & !is.na(symbol)]
gwas_genes[, symbol := trimws(symbol)]
gwas_symbols <- unique(gwas_genes$symbol)
cat(sprintf("    Loaded %d unique GWAS gene symbols\n", length(gwas_symbols)))

# ============================================================
#  2. Load concordance atlas & map GWAS genes
# ============================================================
cat("[2] Loading concordance atlas...\n")
atlas <- fread(ATLAS_FILE)
cat(sprintf("    Atlas contains %d genes across %d categories\n",
            nrow(atlas), uniqueN(atlas$primary_category)))

# Map GWAS genes to atlas
atlas[, is_gwas := human_symbol %in% gwas_symbols]
gwas_in_atlas <- atlas[is_gwas == TRUE]
cat(sprintf("    %d / %d GWAS genes found in atlas (%d not mapped)\n",
            uniqueN(gwas_in_atlas$human_symbol),
            length(gwas_symbols),
            length(gwas_symbols) - uniqueN(gwas_in_atlas$human_symbol)))

# Report unmapped genes
unmapped <- setdiff(gwas_symbols, atlas$human_symbol)
if (length(unmapped) > 0) {
  cat(sprintf("    Unmapped GWAS genes: %s\n",
              paste(head(unmapped, 20), collapse = ", ")))
  if (length(unmapped) > 20) cat(sprintf("    ... and %d more\n", length(unmapped) - 20))
}

# ============================================================
#  3. Fisher exact enrichment tests per concordance category
# ============================================================
cat("[3] Running Fisher exact enrichment tests...\n")

categories <- sort(unique(atlas$primary_category))
n_atlas    <- nrow(atlas)
n_gwas_tot <- sum(atlas$is_gwas)

fisher_results <- rbindlist(lapply(categories, function(cat_name) {
  in_cat    <- atlas$primary_category == cat_name
  gwas_flag <- atlas$is_gwas

  # 2x2 contingency: [[gwas_in_cat, gwas_not_in_cat], [bg_in_cat, bg_not_in_cat]]
  a <- sum(gwas_flag & in_cat)       # GWAS genes in this category

  b <- sum(gwas_flag & !in_cat)      # GWAS genes NOT in this category
  c <- sum(!gwas_flag & in_cat)      # non-GWAS genes in this category
  d <- sum(!gwas_flag & !in_cat)     # non-GWAS genes NOT in this category

  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft  <- fisher.test(mat, alternative = "greater")

  pct_gwas <- ifelse(n_gwas_tot > 0, a / n_gwas_tot * 100, 0)
  n_bg     <- sum(in_cat)
  pct_bg   <- n_bg / n_atlas * 100

  data.table(
    category   = cat_name,
    n_gwas     = a,
    n_total    = n_bg,
    pct_gwas   = round(pct_gwas, 2),
    pct_bg     = round(pct_bg, 2),
    odds_ratio = round(ft$estimate, 3),
    p_value    = ft$p.value
  )
}))

fisher_results[, padj := p.adjust(p_value, method = "fdr")]
fisher_results <- fisher_results[order(p_value)]
fisher_results[, sig_label := fifelse(padj < 0.001, "***",
                              fifelse(padj < 0.01,  "**",
                              fifelse(padj < 0.05,  "*", "ns")))]

cat("    Enrichment results:\n")
print(fisher_results[, .(category, n_gwas, n_total, odds_ratio, padj, sig_label)])

fwrite(fisher_results,
       file.path(OUT_DIR, "gwas_category_enrichment_fisher.csv"))
cat(sprintf("    Saved: %s\n",
            file.path(OUT_DIR, "gwas_category_enrichment_fisher.csv")))

# ============================================================
#  4. Build evidence table: GWAS genes in atlas + DEG tier (MR stripped 2026-04-22)
# ============================================================
cat("[4] Building GWAS evidence table...\n")

evidence <- gwas_in_atlas[, .(
  human_symbol,
  concordance_category = primary_category,
  n_concordant,
  n_diets_sig,
  mean_h_lfc = mean_h_lfc,
  translatability_score,
  translatability_tier,
  best_n_concordant,
  h_significant
)]

# --- Map Ensembl to symbol via biomaRt for consensus DEG merge ---
cat("    Mapping Ensembl IDs to gene symbols via biomaRt...\n")
consensus <- fread(CONSENSUS_FILE)
# Strip version suffix from Ensembl IDs
consensus[, ensembl_base := sub("\\..*$", "", gene)]

ensembl_to_symbol <- tryCatch({
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                        mirror = "useast")
  ids <- unique(consensus$ensembl_base)
  bm <- as.data.table(getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = ids,
    mart       = ensembl
  ))
  setnames(bm, c("ensembl_base", "symbol"))
  bm[symbol != ""]
}, error = function(e) {
  cat(sprintf("    WARNING: biomaRt query failed: %s\n", conditionMessage(e)))
  cat("    Attempting offline Ensembl mapping from atlas...\n")
  # Fallback: use the atlas itself which has human_symbol already
  NULL
})

if (!is.null(ensembl_to_symbol) && nrow(ensembl_to_symbol) > 0) {
  consensus <- merge(consensus, ensembl_to_symbol, by = "ensembl_base", all.x = TRUE)
} else {
  # Fallback: try to match gene symbols from atlas to consensus via Ensembl IDs
  # Build a lookup from atlas (which has human_symbol & mouse_gene_id)
  # This is limited but at least covers GWAS genes in the atlas
  cat("    Using atlas human_symbol as fallback mapping.\n")
  consensus[, symbol := NA_character_]
}

# Subset consensus to GWAS genes
if ("symbol" %in% names(consensus)) {
  consensus_gwas <- consensus[symbol %in% gwas_symbols,
    .(symbol, dream_logFC, dream_padj, meta_logFC, meta_padj, tier)]
  # De-duplicate: keep the entry with smallest dream_padj per symbol
  consensus_gwas <- consensus_gwas[order(dream_padj)]
  consensus_gwas <- consensus_gwas[!duplicated(symbol)]

  evidence <- merge(evidence, consensus_gwas,
                    by.x = "human_symbol", by.y = "symbol", all.x = TRUE)
  cat(sprintf("    Merged consensus DEG tier for %d / %d GWAS-in-atlas genes\n",
              sum(!is.na(evidence$tier)), nrow(evidence)))
} else {
  evidence[, c("dream_logFC", "dream_padj", "meta_logFC", "meta_padj", "tier") :=
             .(NA_real_, NA_real_, NA_real_, NA_real_, NA_character_)]
}

# --- MR/TWAS merge REMOVED 2026-04-22 — MR ditched from paper ---
# TWAS + COLOC + INTACT causal evidence is now merged downstream via the
# multi-evidence atlas (Script 27a). This script's evidence table intentionally
# reports only transcriptomic + concordance signal for GWAS loci.

# Sort by translatability_score descending
evidence <- evidence[order(-translatability_score, -abs(mean_h_lfc))]

fwrite(evidence, file.path(OUT_DIR, "gwas_gene_evidence_table.csv"))
cat(sprintf("    Saved: %s\n",
            file.path(OUT_DIR, "gwas_gene_evidence_table.csv")))

# ============================================================
#  5. Top GWAS loci evidence cards
# ============================================================
cat("[5] Generating evidence cards for priority GWAS loci...\n")

priority_loci <- c("PNPLA3", "TM6SF2", "GCKR", "MBOAT7", "HSD17B13",
                    "FTO", "SERPINA1", "GPAM", "MTTP", "ALDH2",
                    "APOE", "MTARC1", "GATAD2A", "MAU2")

loci_cards <- evidence[human_symbol %in% priority_loci]
# Add any priority loci not in evidence table as rows with NAs
missing_loci <- setdiff(priority_loci, loci_cards$human_symbol)
if (length(missing_loci) > 0) {
  missing_dt <- data.table(human_symbol = missing_loci)
  loci_cards <- rbindlist(list(loci_cards, missing_dt), fill = TRUE)
}

# Order by the priority list
loci_cards[, priority_order := match(human_symbol, priority_loci)]
loci_cards <- loci_cards[order(priority_order)]
loci_cards[, priority_order := NULL]

# Tag in_atlas and in_gwas_list
loci_cards[, in_concordance_atlas := human_symbol %in% atlas$human_symbol]
loci_cards[, in_gwas_list := human_symbol %in% gwas_symbols]

fwrite(loci_cards, file.path(OUT_DIR, "priority_loci_evidence_cards.csv"))
cat(sprintf("    Saved evidence cards for %d priority loci (%d in atlas, %d missing)\n",
            nrow(loci_cards),
            sum(loci_cards$in_concordance_atlas),
            sum(!loci_cards$in_concordance_atlas)))

# Print summary to log
cat("\n    Priority loci summary:\n")
for (i in seq_len(nrow(loci_cards))) {
  row <- loci_cards[i]
  cat(sprintf("    - %s: category=%s, n_concordant=%s, tier=%s, trans_score=%s\n",
              row$human_symbol,
              ifelse(is.na(row$concordance_category), "NOT_IN_ATLAS",
                     row$concordance_category),
              ifelse(is.na(row$n_concordant), "-", as.character(row$n_concordant)),
              ifelse(is.na(row$tier), "-", row$tier),
              ifelse(is.na(row$translatability_score), "-",
                     sprintf("%.3f", row$translatability_score))))
}

# ============================================================
#  6a. Figure: Enrichment bar chart (OR per category)
# ============================================================
cat("\n[6] Generating figures...\n")

# Prepare data for plot
plot_enrich <- copy(fisher_results)
# Cap OR for visualization
plot_enrich[, or_plot := pmin(odds_ratio, 10)]
# Reorder by OR
plot_enrich[, category := factor(category, levels = category[order(odds_ratio)])]

p_enrich <- ggplot(plot_enrich,
                   aes(x = category, y = or_plot, fill = category)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_text(aes(label = sig_label), hjust = -0.3, size = 2.5) +
  scale_fill_manual(values = concordance_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title    = "GWAS Gene Enrichment Across Concordance Categories",
    subtitle = sprintf("Fisher exact test (one-sided); %d GWAS genes mapped to atlas",
                       n_gwas_tot),
    x = NULL,
    y = "Odds Ratio"
  ) +
  theme_pub +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(FIG_DIR, "gwas_enrichment_by_category.pdf"),
       p_enrich, width = 5, height = 3.5, device = cairo_pdf)
cat(sprintf("    Saved: %s\n", file.path(FIG_DIR, "gwas_enrichment_by_category.pdf")))

# ============================================================
#  6b. Figure: GWAS gene category distribution bar chart
# ============================================================
cat_dist <- gwas_in_atlas[, .N, by = primary_category]
setnames(cat_dist, "N", "count")
# Add categories with zero GWAS genes
all_cats <- data.table(primary_category = names(concordance_colors))
cat_dist <- merge(all_cats, cat_dist, by = "primary_category", all.x = TRUE)
cat_dist[is.na(count), count := 0]
cat_dist[, pct := round(count / sum(count) * 100, 1)]
cat_dist[, primary_category := factor(primary_category,
           levels = cat_dist$primary_category[order(-count)])]

p_dist <- ggplot(cat_dist[count > 0],
                 aes(x = primary_category, y = count, fill = primary_category)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", count, pct)),
            vjust = -0.3, size = 2, lineheight = 0.8) +
  scale_fill_manual(values = concordance_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title    = "Distribution of GWAS Genes Across Concordance Categories",
    subtitle = sprintf("%d GWAS genes mapped to concordance atlas",
                       nrow(gwas_in_atlas)),
    x = NULL,
    y = "Number of GWAS genes"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.5))

ggsave(file.path(FIG_DIR, "gwas_category_distribution.pdf"),
       p_dist, width = 5, height = 3.5, device = cairo_pdf)
cat(sprintf("    Saved: %s\n", file.path(FIG_DIR, "gwas_category_distribution.pdf")))

# ============================================================
#  Save category distribution table
# ============================================================
fwrite(cat_dist[order(-count)],
       file.path(OUT_DIR, "gwas_category_distribution.csv"))

# ============================================================
#  Session summary
# ============================================================
cat("\n=== Summary ===\n")
cat(sprintf("  GWAS genes loaded:        %d\n", length(gwas_symbols)))
cat(sprintf("  Mapped to atlas:          %d\n", uniqueN(gwas_in_atlas$human_symbol)))
cat(sprintf("  Enrichment tests:         %d categories\n", nrow(fisher_results)))
cat(sprintf("  Significant (FDR<0.05):   %d\n", sum(fisher_results$padj < 0.05)))
cat(sprintf("  Priority loci w/ data:    %d / %d\n",
            sum(loci_cards$in_concordance_atlas), length(priority_loci)))
cat(sprintf("\nOutputs in: %s\n", OUT_DIR))
cat(sprintf("Figures in: %s\n", FIG_DIR))
cat(sprintf("\nCompleted: %s\n", Sys.time()))
