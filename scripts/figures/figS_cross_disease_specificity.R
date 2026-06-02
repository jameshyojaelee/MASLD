#!/usr/bin/env Rscript
# ===========================================================================
# figS_cross_disease_specificity.R
# Gap #5: Cross-Disease Specificity — Are MASLD DEGs generic liver disease
#         genes or MASLD-specific?
#
# Compares MASLD dream DEGs against ALD, HCV, and generic fibrosis gene sets
# to quantify disease specificity.
#
# Gene set sources:
#   - MASLD DEGs: dream mega-analysis (padj<0.05, |logFC|>0.5)
#   - ALD: Curated from Argemi et al. 2019 Hepatology + msigdbr CGP
#   - HCV: Curated ISG signature + msigdbr WP_HEPATITIS_C pathway
#   - Generic fibrosis: Govaere 25-gene panel + KEGG/Reactome ECM/collagen
#   - NAFLD pathway: msigdbr WP_NONALCOHOLIC_FATTY_LIVER_DISEASE (155 genes)
#
# Output:
#   figures/supplementary/figS_sensitivity/figS_cross_disease_specificity.pdf
#   figures/supplementary/figS_sensitivity/cross_disease_specificity_data.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(grid)
  library(msigdbr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUTDIR <- FIGS_SENS_DIR
dir.create(file.path(OUTDIR, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("=== Cross-Disease Specificity Analysis ===\n")

# ============================================================================
# 1. Load MASLD DEGs
# ============================================================================
dream_file <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                        "results/integration/dream_results_ashr.csv")
dream <- fread(dream_file)
cat("Loaded dream results:", nrow(dream), "genes\n")

# Define MASLD DEGs
masld_degs <- dream[padj < 0.05 & abs(logFC) > 0.5, symbol]
masld_degs <- masld_degs[!is.na(masld_degs) & masld_degs != ""]
masld_degs <- unique(masld_degs)
cat("MASLD DEGs (padj<0.05, |logFC|>0.5):", length(masld_degs), "\n")

# Also split by direction
masld_up <- dream[padj < 0.05 & logFC > 0.5, symbol]
masld_up <- unique(masld_up[!is.na(masld_up) & masld_up != ""])
masld_dn <- dream[padj < 0.05 & logFC < -0.5, symbol]
masld_dn <- unique(masld_dn[!is.na(masld_dn) & masld_dn != ""])

# Background: all tested genes
all_tested <- dream[!is.na(symbol) & symbol != "", unique(symbol)]
cat("Background (all tested):", length(all_tested), "\n\n")

# ============================================================================
# 2. Assemble disease gene sets
# ============================================================================

# --- 2a. ALD signature ---
# Core ALD genes from Argemi et al. 2019 Hepatology (alcoholic hepatitis
# transcriptomics, 77 patients) + Affoh et al. 2020 + literature consensus.
# Includes: alcohol metabolism, neutrophil chemokines, acute phase, ER stress
ald_curated <- c(
  # Alcohol metabolism (downregulated in ALD)
  "ADH1B", "ADH1C", "ADH4", "ADH6", "ALDH2", "CYP2E1", "CYP3A4", "CYP2A6",
  "CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP2B6", "HSD17B6", "HSD17B13",
  # Neutrophil chemokines / inflammation (upregulated in ALD)
  "CXCL1", "CXCL5", "CXCL6", "CXCL8", "CCL2", "CCL20", "IL8",
  "TNF", "IL6", "IL1B", "IL1A", "TNFRSF12A",
  # Acute phase / liver damage markers
  "SAA1", "SAA2", "LCN2", "LBP", "CRP", "HP", "ORM1", "ORM2",
  "SERPINA3", "A2M",
  # Lipid metabolism (dysregulated in ALD)
  "AKR1B10", "FABP4", "FASN", "SCD", "ACLY",
  # Liver injury / regeneration
  "SPP1", "KRT23", "KRT7", "ANXA2", "ANXA1", "MMP7", "MMP9",
  # Fibrosis / ECM (shared with other liver diseases)
  "COL1A1", "COL1A2", "COL3A1", "ACTA2", "TIMP1", "TGFB1",
  "LOX", "LOXL2", "THY1",
  # ER stress
  "DDIT3", "ATF4", "XBP1", "HSPA5", "EIF2AK3",
  # Oxidative stress
  "NQO1", "GCLC", "GCLM", "SOD2", "GPX1",
  # Complement (activated in ALD)
  "C3", "C4A", "C4B", "C5", "CFB", "CFH",
  # Interferon response (modest in ALD unlike HCV)
  "ISG15", "IFIT1", "MX1"
)

# --- 2b. HCV signature ---
# HCV liver transcriptome is dominated by interferon-stimulated genes (ISGs).
# Core ISG signature from Hoshida 2013, Liang 2013, Sarasin-Filipowicz 2008,
# Bolen 2014, and msigdbr WP pathway.
hcv_curated <- c(
  # Core ISG signature (hallmark of HCV liver)
  "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFIT5",
  "IFI44", "IFI44L", "IFI6", "IFI27", "IFI35",
  "MX1", "MX2",
  "OAS1", "OAS2", "OAS3", "OASL",
  "RSAD2", "HERC5", "HERC6",
  "DDX58", "DDX60", "DDX60L",
  "STAT1", "STAT2", "IRF7", "IRF9", "IRF1",
  "CXCL10", "CXCL9", "CXCL11",
  # Antiviral effectors
  "GBP1", "GBP2", "GBP4", "GBP5",
  "TRIM22", "TRIM25", "TRIM5",
  "PLSCR1", "EIF2AK2", "ADAR",
  "USP18", "BST2", "LGALS3BP",
  # T-cell / immune infiltrate (chronic HCV)
  "CD8A", "CD8B", "GZMB", "PRF1", "GNLY",
  "CXCR3", "CCL5", "LAG3", "PDCD1", "HAVCR2",
  # HCV-specific fibrosis markers
  "COL1A1", "COL1A2", "COL3A1", "ACTA2", "TIMP1",
  "LOX", "LOXL2", "MMP2",
  # Lipid metabolism (HCV hijacks lipid pathways)
  "LDLR", "HMGCR", "SREBF2", "PCSK9",
  # Apoptosis / liver damage
  "TRAIL", "FAS", "CASP3", "BCL2"
)

# Supplement with msigdbr WP_HEPATITIS_C_AND_HEPATOCELLULAR_CARCINOMA
msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
wp_hcv <- msig_c2[msig_c2$gs_name == "WP_HEPATITIS_C_AND_HEPATOCELLULAR_CARCINOMA",
                   "gene_symbol", drop = FALSE]
hcv_curated <- unique(c(hcv_curated, wp_hcv$gene_symbol))

# --- 2c. Generic fibrosis gene set ---
# Genes involved in fibrosis across ALL chronic liver diseases (not disease-specific)
# Sources: Govaere 2025 (25-gene panel), KEGG ECM, Reactome collagen,
# and curated core fibrosis genes from Tsuchida & Friedman 2017
govaere_25 <- c("AKR1B10", "DUSP6", "GDF15", "THBS2", "A2M", "CDH2",
                "COL1A1", "COL3A1", "COL4A1", "COL4A2", "COL6A3", "DCN",
                "FBN1", "FSTL1", "IGFBP7", "LUM", "MFAP4", "MMP2",
                "POSTN", "SPARC", "SPP1", "TAGLN", "THY1", "TIMP1", "VCAN")

fibrosis_core <- c(
  # Collagens (pan-fibrosis)
  "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2",
  "COL6A1", "COL6A2", "COL6A3", "COL14A1", "COL15A1", "COL18A1",
  # ECM remodeling
  "FN1", "FBN1", "ELN", "POSTN", "SPARC", "THBS1", "THBS2",
  "TNC", "VCAN", "DCN", "BGN", "LUM", "FMOD",
  # Matrix metalloproteinases & inhibitors
  "MMP2", "MMP7", "MMP9", "MMP14", "TIMP1", "TIMP2", "TIMP3",
  # Crosslinking
  "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4",
  # TGF-beta pathway (master fibrosis driver)
  "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
  "SMAD2", "SMAD3", "SMAD4", "SMAD7",
  # Stellate cell activation
  "ACTA2", "PDGFRB", "PDGFRA", "DES", "VIM",
  "CTGF", "PDGFB", "PDGFA",
  # Integrins
  "ITGB1", "ITGAV", "ITGB6",
  # HSC-activation markers
  "THY1", "MFAP4", "FSTL1", "IGFBP7", "TAGLN",
  "SERPINE1", "CDH2", "LGALS1",
  # Inflammatory fibrosis
  "CCL2", "CCR2", "IL13", "IL4", "IL33",
  # Govaere markers (overlap)
  "A2M", "GDF15", "AKR1B10", "DUSP6", "SPP1"
)
fibrosis_genes <- unique(c(govaere_25, fibrosis_core))

# --- 2d. NAFLD pathway (msigdbr) ---
wp_nafld <- msig_c2[msig_c2$gs_name == "WP_NONALCOHOLIC_FATTY_LIVER_DISEASE",
                    "gene_symbol", drop = FALSE]
nafld_pathway <- unique(wp_nafld$gene_symbol)

# Restrict all sets to genes tested in our data
ald_genes   <- intersect(ald_curated, all_tested)
hcv_genes   <- intersect(hcv_curated, all_tested)
fib_genes   <- intersect(fibrosis_genes, all_tested)
nafld_genes <- intersect(nafld_pathway, all_tested)

cat("Gene set sizes (after restricting to tested genes):\n")
cat("  ALD:", length(ald_genes), "/", length(ald_curated), "\n")
cat("  HCV:", length(hcv_genes), "/", length(hcv_curated), "\n")
cat("  Generic fibrosis:", length(fib_genes), "/", length(fibrosis_genes), "\n")
cat("  WP NAFLD pathway:", length(nafld_genes), "/", length(nafld_pathway), "\n\n")

# ============================================================================
# 3. Overlap and specificity analysis
# ============================================================================

# Intersections
masld_and_ald <- intersect(masld_degs, ald_genes)
masld_and_hcv <- intersect(masld_degs, hcv_genes)
masld_and_fib <- intersect(masld_degs, fib_genes)
ald_and_hcv   <- intersect(ald_genes, hcv_genes)
shared_all3   <- Reduce(intersect, list(masld_degs, ald_genes, hcv_genes))

# MASLD-specific: in our DEGs but NOT in ALD or HCV
masld_specific <- setdiff(masld_degs, union(ald_genes, hcv_genes))
masld_not_fib  <- setdiff(masld_degs, fib_genes)

cat("Overlap counts:\n")
cat("  MASLD ∩ ALD:", length(masld_and_ald), "\n")
cat("  MASLD ∩ HCV:", length(masld_and_hcv), "\n")
cat("  MASLD ∩ Fibrosis:", length(masld_and_fib), "\n")
cat("  ALD ∩ HCV:", length(ald_and_hcv), "\n")
cat("  All three:", length(shared_all3), "\n")
cat("  MASLD-specific (not ALD/HCV):", length(masld_specific),
    sprintf("(%.1f%%)\n", 100 * length(masld_specific) / length(masld_degs)))
cat("  MASLD not in fibrosis:", length(masld_not_fib),
    sprintf("(%.1f%%)\n\n", 100 * length(masld_not_fib) / length(masld_degs)))

# Jaccard indices
jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

jac_masld_ald <- jaccard(masld_degs, ald_genes)
jac_masld_hcv <- jaccard(masld_degs, hcv_genes)
jac_masld_fib <- jaccard(masld_degs, fib_genes)
jac_ald_hcv   <- jaccard(ald_genes, hcv_genes)
jac_ald_fib   <- jaccard(ald_genes, fib_genes)
jac_hcv_fib   <- jaccard(hcv_genes, fib_genes)

cat("Jaccard indices:\n")
cat("  MASLD vs ALD:", round(jac_masld_ald, 4), "\n")
cat("  MASLD vs HCV:", round(jac_masld_hcv, 4), "\n")
cat("  MASLD vs Fibrosis:", round(jac_masld_fib, 4), "\n")
cat("  ALD vs HCV:", round(jac_ald_hcv, 4), "\n")
cat("  ALD vs Fibrosis:", round(jac_ald_fib, 4), "\n")
cat("  HCV vs Fibrosis:", round(jac_hcv_fib, 4), "\n\n")

# Fisher exact tests (enrichment of MASLD DEGs in each disease set)
fisher_test <- function(degs, disease_set, background) {
  a <- length(intersect(degs, disease_set))
  b <- length(setdiff(degs, disease_set))
  c <- length(setdiff(disease_set, degs))
  d <- length(setdiff(background, union(degs, disease_set)))
  mat <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  list(overlap = a, OR = ft$estimate, p = ft$p.value,
       ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2])
}

ft_ald <- fisher_test(masld_degs, ald_genes, all_tested)
ft_hcv <- fisher_test(masld_degs, hcv_genes, all_tested)
ft_fib <- fisher_test(masld_degs, fib_genes, all_tested)

cat("Fisher enrichment (MASLD DEGs in disease sets):\n")
cat(sprintf("  ALD: OR=%.2f, p=%.2e, overlap=%d\n", ft_ald$OR, ft_ald$p, ft_ald$overlap))
cat(sprintf("  HCV: OR=%.2f, p=%.2e, overlap=%d\n", ft_hcv$OR, ft_hcv$p, ft_hcv$overlap))
cat(sprintf("  Fibrosis: OR=%.2f, p=%.2e, overlap=%d\n\n", ft_fib$OR, ft_fib$p, ft_fib$overlap))

# ============================================================================
# 4. Classify MASLD DEGs by disease specificity
# ============================================================================
classify_gene <- function(g) {
  in_ald <- g %in% ald_genes
  in_hcv <- g %in% hcv_genes
  in_fib <- g %in% fib_genes
  if (in_ald & in_hcv & in_fib) return("Shared (all)")
  if (in_ald & in_hcv) return("ALD+HCV shared")
  if (in_ald & in_fib) return("ALD+Fibrosis")
  if (in_hcv & in_fib) return("HCV+Fibrosis")
  if (in_ald) return("ALD only")
  if (in_hcv) return("HCV only")
  if (in_fib) return("Fibrosis only")
  return("MASLD-specific")
}

masld_class <- data.table(
  gene = masld_degs,
  logFC = dream[match(masld_degs, symbol), logFC],
  padj  = dream[match(masld_degs, symbol), padj],
  category = vapply(masld_degs, classify_gene, character(1))
)
masld_class[, direction := ifelse(logFC > 0, "Up", "Down")]

cat("MASLD DEG classification:\n")
print(masld_class[, .N, by = category][order(-N)])

# ============================================================================
# 5. Build output data
# ============================================================================
# Per-gene detail
out_dt <- copy(masld_class)
out_dt[, in_ALD := gene %in% ald_genes]
out_dt[, in_HCV := gene %in% hcv_genes]
out_dt[, in_fibrosis := gene %in% fib_genes]
out_dt[, in_WP_NAFLD := gene %in% nafld_genes]

fwrite(out_dt, file.path(OUTDIR, "cross_disease_specificity_data.csv"))
cat("\nSaved per-gene data:", nrow(out_dt), "rows\n")

# Summary statistics
summary_dt <- data.table(
  comparison = c("MASLD_vs_ALD", "MASLD_vs_HCV", "MASLD_vs_Fibrosis",
                 "ALD_vs_HCV", "ALD_vs_Fibrosis", "HCV_vs_Fibrosis"),
  jaccard = c(jac_masld_ald, jac_masld_hcv, jac_masld_fib,
              jac_ald_hcv, jac_ald_fib, jac_hcv_fib),
  overlap = c(length(masld_and_ald), length(masld_and_hcv), length(masld_and_fib),
              length(ald_and_hcv), length(intersect(ald_genes, fib_genes)),
              length(intersect(hcv_genes, fib_genes))),
  set1_size = c(rep(length(masld_degs), 3), length(ald_genes),
                length(ald_genes), length(hcv_genes)),
  set2_size = c(length(ald_genes), length(hcv_genes), length(fib_genes),
                length(hcv_genes), length(fib_genes), length(fib_genes))
)
fwrite(summary_dt, file.path(OUTDIR, "cross_disease_jaccard_summary.csv"))

# ============================================================================
# 6. FIGURE PANELS
# ============================================================================
cat("\n=== Generating figure panels ===\n")

# --- Panel (a): UpSet-style visualization using ComplexHeatmap ---
# Build a binary membership matrix for all genes in any set
all_union <- unique(c(masld_degs, ald_genes, hcv_genes, fib_genes))
m <- data.frame(
  MASLD = as.integer(all_union %in% masld_degs),
  ALD   = as.integer(all_union %in% ald_genes),
  HCV   = as.integer(all_union %in% hcv_genes),
  Fibrosis = as.integer(all_union %in% fib_genes),
  row.names = all_union
)
comb_mat <- make_comb_mat(m, mode = "distinct")

# Filter to combinations with at least 3 members for readability
comb_mat_filtered <- comb_mat[comb_size(comb_mat) >= 3]

upset_colors <- c(
  "MASLD"    = masld_colors$masld,
  "ALD"      = "#FF8F00",
  "HCV"      = "#2E7D32",
  "Fibrosis" = "#5D4037"
)

# --- Panel (b): Jaccard similarity heatmap ---
sets <- list(MASLD = masld_degs, ALD = ald_genes, HCV = hcv_genes,
             Fibrosis = fib_genes)
n <- length(sets)
jac_mat <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
pval_mat <- matrix(1, n, n, dimnames = list(names(sets), names(sets)))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i == j) {
      jac_mat[i, j] <- 1
      pval_mat[i, j] <- 0
    } else {
      jac_mat[i, j] <- jaccard(sets[[i]], sets[[j]])
      ft <- fisher_test(sets[[i]], sets[[j]], all_tested)
      pval_mat[i, j] <- ft$p
    }
  }
}

# Format p-value labels
pval_labels <- matrix("", n, n)
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i == j) {
      pval_labels[i, j] <- "1.000"
    } else {
      jv <- sprintf("%.3f", jac_mat[i, j])
      pv <- pval_mat[i, j]
      pstar <- if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else ""
      pval_labels[i, j] <- paste0(jv, pstar)
    }
  }
}

# ggplot heatmap for panel b
jac_long <- data.table(
  row = rep(names(sets), each = n),
  col = rep(names(sets), n),
  jaccard = as.vector(jac_mat),
  label = as.vector(pval_labels)
)
jac_long[, row := factor(row, levels = names(sets))]
jac_long[, col := factor(col, levels = rev(names(sets)))]

pb <- ggplot(jac_long, aes(x = row, y = col, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.2, color = "black") +
  scale_fill_gradient2(low = "white", mid = "#F48FB1", high = "#880E4F",
                       midpoint = 0.15, limits = c(0, 1),
                       name = "Jaccard\nindex") +
  labs(x = NULL, y = NULL, title = "Pairwise Jaccard similarity") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# --- Panel (c): MASLD-specific genes — top by effect size ---
# Get top 30 MASLD-specific genes by |logFC|
top_specific <- masld_class[category == "MASLD-specific"][order(-abs(logFC))][1:min(30, .N)]

pc <- ggplot(top_specific, aes(x = reorder(gene, logFC), y = logFC,
                                fill = direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Up" = masld_colors$up, "Down" = masld_colors$down),
                    name = "Direction") +
  labs(x = NULL, y = "logFC (dream)",
       title = "Top 30 MASLD-specific DEGs",
       subtitle = "Not in ALD, HCV, or fibrosis sets") +
  theme_masld() +
  theme(legend.position = "bottom")

# --- Panel (d): Composition pie/bar ---
cat_summary <- masld_class[, .N, by = category][order(-N)]
cat_summary[, pct := round(100 * N / sum(N), 1)]

# Simplify categories for display
cat_summary[, display := category]
# Merge small categories
small_cats <- cat_summary[N < 10, category]
if (length(small_cats) > 0) {
  cat_summary[category %in% small_cats, display := "Other shared"]
  cat_summary <- cat_summary[, .(N = sum(N)), by = display]
  cat_summary[, pct := round(100 * N / sum(N), 1)]
}

# Order: MASLD-specific first, then by size
cat_summary[, display := factor(display,
  levels = c("MASLD-specific",
             sort(setdiff(cat_summary$display, "MASLD-specific"), decreasing = TRUE)))]

# Disease-specificity palette
spec_colors <- c(
  "MASLD-specific" = masld_colors$masld,
  "Fibrosis only"  = "#5D4037",
  "ALD only"       = "#FF8F00",
  "HCV only"       = "#2E7D32",
  "Shared (all)"   = "#9E9E9E",
  "ALD+Fibrosis"   = "#A1887F",
  "HCV+Fibrosis"   = "#66BB6A",
  "ALD+HCV shared" = "#FFC107",
  "Other shared"   = "#BDBDBD"
)

pd <- ggplot(cat_summary, aes(x = display, y = N, fill = display)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", N, pct)),
            vjust = -0.3, size = 2) +
  scale_fill_manual(values = spec_colors, guide = "none") +
  labs(x = NULL, y = "Number of MASLD DEGs",
       title = "Disease specificity of MASLD DEGs") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# ============================================================================
# 7. Compose and save figure
# ============================================================================

# Panel (a) is ComplexHeatmap (grid), panels b-d are ggplot
# Strategy: save UpSet as separate grob, then combine with patchwork via wrap_elements

# Save UpSet to temporary PDF then use as grob
upset_grob <- grid.grabExpr({
  draw(UpSet(
    comb_mat_filtered,
    set_order = c("MASLD", "ALD", "HCV", "Fibrosis"),
    comb_order = order(-comb_size(comb_mat_filtered)),
    top_annotation = upset_top_annotation(comb_mat_filtered,
      height = unit(3, "cm"),
      annotation_name_side = "left",
      annotation_name_gp = gpar(fontsize = 6),
      axis_param = list(gp = gpar(fontsize = 5))),
    left_annotation = upset_left_annotation(comb_mat_filtered,
      width = unit(2, "cm"),
      annotation_name_gp = gpar(fontsize = 6),
      axis_param = list(gp = gpar(fontsize = 5))),
    row_names_gp = gpar(fontsize = 7),
    column_title = "Gene set intersections",
    column_title_gp = gpar(fontsize = 8, fontface = "bold")
  ))
}, width = 5, height = 4)

pa_wrapped <- wrap_elements(upset_grob)

# Compose: top row = UpSet + Jaccard, bottom row = specific genes + composition
fig <- (pa_wrapped | pb) /
       (pc | pd) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

out_pdf <- file.path(OUTDIR, "figS_cross_disease_specificity.pdf")
save_fig(fig, out_pdf, width = fig_full_width, height = 8)
cat("\nSaved figure:", out_pdf, "\n")

# ============================================================================
# 8. Print final summary
# ============================================================================
cat("\n=== SUMMARY ===\n")
cat(sprintf("Total MASLD DEGs: %d\n", length(masld_degs)))
cat(sprintf("MASLD-specific (not in ALD/HCV/Fibrosis): %d (%.1f%%)\n",
            sum(masld_class$category == "MASLD-specific"),
            100 * sum(masld_class$category == "MASLD-specific") / length(masld_degs)))
cat(sprintf("Shared with ALD: %d (%.1f%%) | Jaccard=%.3f | Fisher OR=%.2f, p=%.1e\n",
            ft_ald$overlap, 100 * ft_ald$overlap / length(masld_degs),
            jac_masld_ald, ft_ald$OR, ft_ald$p))
cat(sprintf("Shared with HCV: %d (%.1f%%) | Jaccard=%.3f | Fisher OR=%.2f, p=%.1e\n",
            ft_hcv$overlap, 100 * ft_hcv$overlap / length(masld_degs),
            jac_masld_hcv, ft_hcv$OR, ft_hcv$p))
cat(sprintf("Shared with Fibrosis: %d (%.1f%%) | Jaccard=%.3f | Fisher OR=%.2f, p=%.1e\n",
            ft_fib$overlap, 100 * ft_fib$overlap / length(masld_degs),
            jac_masld_fib, ft_fib$OR, ft_fib$p))
cat(sprintf("Shared across all 3 diseases: %d genes\n", length(shared_all3)))
if (length(shared_all3) > 0 & length(shared_all3) <= 30) {
  cat("  Genes:", paste(sort(shared_all3), collapse = ", "), "\n")
}

cat("\nGene set sources:\n")
cat("  ALD: Curated (Argemi 2019 + literature), N=", length(ald_genes), "\n")
cat("  HCV: Curated (ISG core + msigdbr WP), N=", length(hcv_genes), "\n")
cat("  Fibrosis: Govaere 25 + curated ECM/TGFb, N=", length(fib_genes), "\n")

cat("\nDone.\n")
