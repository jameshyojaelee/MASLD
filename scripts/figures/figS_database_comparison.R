##############################################################################
# Supplementary Figure: External Database Comparison
# 4 panels: (a) UpSet overlap, (b) Recovery bar chart,
#           (c) Atlas-unique top genes, (d) Statistics summary text
#
# Databases:
#   1. KEGG hsa04932 NAFLD pathway (157 genes, via KEGGREST)
#   2. WikiPathways NAFLD (155 genes, via msigdbr)
#   3. Open Targets NAFLD/FLD (5,946 genes, pre-downloaded)
#   4. Atlas DEGs (lfsr<0.05, |shrunk_logFC|>0.3)
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
  library(ComplexHeatmap)
  library(grid)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_SENS_DIR
dir.create(file.path(OUT_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("=== Loading data ===\n")

# ---- 1. Atlas (multi-evidence) ----
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
cat("Atlas genes:", nrow(atlas), "\n")

# DEGs: C2 canonical Tier-1 (bulk_lfsr<0.05, |bulk_shrunk_logFC|>0.3)
atlas[, is_deg := is_dream_deg(atlas)]
deg_genes <- atlas[is_deg == TRUE, human_symbol]
cat("Atlas DEGs:", length(deg_genes), "\n")

# COLOC support: any coloc PP4 > 0.5
coloc_cols <- grep("coloc_pp4", names(atlas), value = TRUE)
if (length(coloc_cols) > 0) {
  atlas[, has_coloc := apply(.SD, 1, function(x) any(!is.na(x) & x > 0.5)), .SDcols = coloc_cols]
} else {
  atlas[, has_coloc := FALSE]
}
cat("Genes with COLOC PP4>0.5:", sum(atlas$has_coloc, na.rm = TRUE), "\n")

# ---- 2. KEGG hsa04932 NAFLD pathway ----
kegg_genes <- tryCatch({
  library(KEGGREST)
  pathway <- keggGet("hsa04932")
  raw <- pathway[[1]]$GENE
  symbols <- raw[seq(2, length(raw), 2)]
  symbols <- sub(";.*", "", symbols)
  # Remove KO entries
  symbols <- symbols[!grepl("^\\[KO:", symbols)]
  symbols
}, error = function(e) {
  cat("KEGGREST failed, using fallback\n")
  # Fallback: well-known KEGG NAFLD genes
  c("PNPLA3", "ADIPOQ", "PPARA", "PPARG", "TNF", "IL6", "IL1B",
    "NFKB1", "TGFB1", "CASP3", "BAX", "BCL2", "INS", "IRS1", "IRS2",
    "AKT1", "AKT2", "PIK3CA", "RXRA", "PPARGC1A", "SIRT1", "MLXIPL",
    "SREBF1", "FASN", "ACACA", "SCD", "DGAT1", "DGAT2", "MTTP", "CYP2E1")
})
cat("KEGG NAFLD genes:", length(kegg_genes), "\n")

# ---- 3. WikiPathways NAFLD ----
wp_genes <- tryCatch({
  library(msigdbr)
  x <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
  nafld <- x[x$gs_name == "WP_NONALCOHOLIC_FATTY_LIVER_DISEASE", ]
  unique(nafld$gene_symbol)
}, error = function(e) {
  cat("msigdbr failed for WP, skipping\n")
  character(0)
})
cat("WikiPathways NAFLD genes:", length(wp_genes), "\n")

# ---- 4. Open Targets NAFLD ----
ot_file <- file.path(BASE, "data/opentargets_nafld_all.csv")
if (file.exists(ot_file)) {
  ot <- fread(ot_file)
  ot_genes <- unique(ot$gene_symbol)
  ot_scores <- ot  # keep scores for later
} else {
  cat("Open Targets file not found, skipping\n")
  ot_genes <- character(0)
  ot_scores <- data.table(gene_symbol = character(0), ot_score = numeric(0))
}
cat("Open Targets NAFLD genes:", length(ot_genes), "\n")

# ---- Build membership matrix ----
all_genes <- unique(c(deg_genes, kegg_genes, wp_genes, ot_genes))
cat("Total unique genes across all sources:", length(all_genes), "\n")

membership <- data.table(
  gene = all_genes,
  `Atlas DEGs` = all_genes %in% deg_genes,
  `KEGG NAFLD` = all_genes %in% kegg_genes,
  `WikiPathways NAFLD` = all_genes %in% wp_genes,
  `Open Targets` = all_genes %in% ot_genes
)

# ---- Panel (a): UpSet plot ----
cat("\n=== Panel (a): UpSet plot ===\n")

# Use ComplexHeatmap UpSet
mat <- as.matrix(membership[, -1])
rownames(mat) <- membership$gene
# Convert logical to numeric
mat_num <- mat * 1L

# Create combination matrix
comb_mat <- make_comb_mat(mat_num, mode = "intersect")

# Generate UpSet as PDF panel
upset_pdf <- file.path(OUT_DIR, "panels", "panel_a_upset.pdf")
pdf(upset_pdf, width = 7, height = 4.5)
upset_plot <- UpSet(
  comb_mat,
  set_order = c("Atlas DEGs", "Open Targets", "KEGG NAFLD", "WikiPathways NAFLD"),
  comb_order = order(comb_size(comb_mat), decreasing = TRUE),
  top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE,
                                         numbers_gp = gpar(fontsize = 6),
                                         height = unit(4, "cm")),
  right_annotation = upset_right_annotation(comb_mat, add_numbers = TRUE,
                                             numbers_gp = gpar(fontsize = 6),
                                             width = unit(3, "cm")),
  row_names_gp = gpar(fontsize = 6)
)
draw(upset_plot)
dev.off()
cat("Saved:", upset_pdf, "\n")
message("[caption] Gene Set Overlaps: Atlas DEGs vs External Databases")

# ---- Panel (b): Recovery bar chart ----
cat("\n=== Panel (b): Recovery bar chart ===\n")

db_list <- list(
  "KEGG NAFLD\n(hsa04932)" = kegg_genes,
  "WikiPathways\nNAFLD" = wp_genes,
  "Open Targets\n(score>0)" = ot_genes,
  "Open Targets\n(score>0.1)" = ot_scores[ot_score > 0.1, gene_symbol]
)

# For each DB: total, in atlas, is DEG, has COLOC
recovery_dt <- rbindlist(lapply(names(db_list), function(db_name) {
  genes <- db_list[[db_name]]
  n_total <- length(genes)
  in_atlas <- genes[genes %in% atlas$human_symbol]
  is_deg <- genes[genes %in% deg_genes]
  has_col <- genes[genes %in% atlas[has_coloc == TRUE, human_symbol]]

  data.table(
    Database = db_name,
    Category = c("Total in DB", "In Atlas", "Atlas DEG", "COLOC PP4>0.5"),
    Count = c(n_total, length(in_atlas), length(is_deg), length(has_col)),
    Pct = c(100,
            round(100 * length(in_atlas) / n_total, 1),
            round(100 * length(is_deg) / n_total, 1),
            round(100 * length(has_col) / n_total, 1))
  )
}))

recovery_dt[, Category := factor(Category,
                                  levels = c("Total in DB", "In Atlas", "Atlas DEG", "COLOC PP4>0.5"))]

bar_colors <- c(
  "Total in DB" = "#BDBDBD",
  "In Atlas" = "#42A5F5",
  "Atlas DEG" = "#C2185B",
  "COLOC PP4>0.5" = "#00695C"
)

p_b <- ggplot(recovery_dt, aes(x = Database, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0(Pct, "%")),
            position = position_dodge(width = 0.8), vjust = -0.3, size = GEOM_TEXT_6PT) +
  scale_fill_manual(values = bar_colors) +
  scale_y_log10(labels = comma, breaks = c(1, 10, 100, 1000, 10000)) +
  labs(x = NULL, y = "Gene count (log scale)") +
  theme_masld() +
  theme(
    axis.text.x = element_text(size = 6, lineheight = 0.9),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 6)
  )
message("[caption] Atlas Recovery of External Database Genes")

# ---- Panel (c): Atlas-unique top genes ----
cat("\n=== Panel (c): Atlas-unique genes ===\n")

# Genes that are DEGs but NOT in any external DB
any_external <- unique(c(kegg_genes, wp_genes, ot_genes))
atlas_unique <- atlas[is_deg == TRUE & !(human_symbol %in% any_external)]
cat("Atlas-unique DEGs (not in any external DB):", nrow(atlas_unique), "\n")

# Top 20 by |logFC|
top_unique <- atlas_unique[order(-abs(bulk_logFC))][1:min(20, nrow(atlas_unique))]
top_unique[, direction := ifelse(bulk_logFC > 0, "Upregulated", "Downregulated")]

p_c <- ggplot(top_unique, aes(x = reorder(human_symbol, abs(bulk_logFC)),
                               y = bulk_logFC, fill = direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = masld_colors$up,
                                "Downregulated" = masld_colors$down)) +
  labs(x = NULL, y = "Dream logFC") +
  theme_masld() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
message(sprintf("[caption] Top 20 Atlas-Unique DEGs (%d DEGs not in any surveyed database)",
                nrow(atlas_unique)))

# ---- Panel (d): Statistics text ----
cat("\n=== Panel (d): Key statistics ===\n")

# Compute statistics
stats_lines <- c()

# KEGG recovery
kegg_in_deg <- sum(kegg_genes %in% deg_genes)
stats_lines <- c(stats_lines,
                 sprintf("KEGG NAFLD (hsa04932): %d/%d genes (%.1f%%) are atlas DEGs",
                         kegg_in_deg, length(kegg_genes),
                         100 * kegg_in_deg / length(kegg_genes)))

# WP recovery
wp_in_deg <- sum(wp_genes %in% deg_genes)
stats_lines <- c(stats_lines,
                 sprintf("WikiPathways NAFLD: %d/%d genes (%.1f%%) are atlas DEGs",
                         wp_in_deg, length(wp_genes),
                         100 * wp_in_deg / length(wp_genes)))

# Open Targets recovery
ot_in_deg <- sum(ot_genes %in% deg_genes)
stats_lines <- c(stats_lines,
                 sprintf("Open Targets NAFLD: %d/%d genes (%.1f%%) are atlas DEGs",
                         ot_in_deg, length(ot_genes),
                         100 * ot_in_deg / length(ot_genes)))

# Atlas-unique
stats_lines <- c(stats_lines,
                 sprintf("Atlas-unique DEGs (not in any DB): %d/%d (%.1f%%)",
                         nrow(atlas_unique), length(deg_genes),
                         100 * nrow(atlas_unique) / length(deg_genes)))

# Union of all external DBs
all_external_genes <- unique(c(kegg_genes, wp_genes, ot_genes))
external_in_atlas <- sum(all_external_genes %in% atlas$human_symbol)
external_are_deg <- sum(all_external_genes %in% deg_genes)
stats_lines <- c(stats_lines,
                 sprintf("Union of all external DBs: %d unique genes", length(all_external_genes)),
                 sprintf("  - %d (%.1f%%) present in atlas", external_in_atlas,
                         100 * external_in_atlas / length(all_external_genes)),
                 sprintf("  - %d (%.1f%%) are atlas DEGs", external_are_deg,
                         100 * external_are_deg / length(all_external_genes)))

# Enrichment: are external DB genes enriched among DEGs?
# Fisher's exact test
in_db_deg <- sum(atlas$human_symbol %in% all_external_genes & atlas$is_deg)
in_db_notdeg <- sum(atlas$human_symbol %in% all_external_genes & !atlas$is_deg)
notin_db_deg <- sum(!(atlas$human_symbol %in% all_external_genes) & atlas$is_deg)
notin_db_notdeg <- sum(!(atlas$human_symbol %in% all_external_genes) & !atlas$is_deg)
fisher_res <- fisher.test(matrix(c(in_db_deg, in_db_notdeg, notin_db_deg, notin_db_notdeg), nrow = 2))
stats_lines <- c(stats_lines,
                 sprintf("Fisher enrichment (DB genes among DEGs): OR=%.2f, P=%.2e",
                         fisher_res$estimate, fisher_res$p.value))

# COLOC enrichment among external DB genes
in_db_coloc <- sum(atlas$human_symbol %in% all_external_genes & atlas$has_coloc)
in_db_nocoloc <- sum(atlas$human_symbol %in% all_external_genes & !atlas$has_coloc)
notin_db_coloc <- sum(!(atlas$human_symbol %in% all_external_genes) & atlas$has_coloc)
notin_db_nocoloc <- sum(!(atlas$human_symbol %in% all_external_genes) & !atlas$has_coloc)
fisher_coloc <- fisher.test(matrix(c(in_db_coloc, in_db_nocoloc, notin_db_coloc, notin_db_nocoloc), nrow = 2))
stats_lines <- c(stats_lines,
                 sprintf("Fisher enrichment (DB genes with COLOC): OR=%.2f, P=%.2e",
                         fisher_coloc$estimate, fisher_coloc$p.value))

# Print all stats
cat("\n--- KEY STATISTICS ---\n")
for (s in stats_lines) cat(s, "\n")
cat("--- END ---\n\n")

# Create text grob for panel (d)
stats_text <- paste(stats_lines, collapse = "\n")
p_d <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = stats_text,
           hjust = 0.5, vjust = 0.5, size = GEOM_TEXT_6PT, family = "mono", lineheight = 1.3) +
  theme_void()
message("[caption] External Database Comparison: Key Statistics")

# ---- Composite figure (panels b-d; panel a is ComplexHeatmap PDF) ----
cat("=== Assembling composite figure ===\n")

composite_pdf <- file.path(OUT_DIR, "figS_database_comparison.pdf")
pdf(composite_pdf, width = 14, height = 14)

# Page 1: UpSet (ComplexHeatmap)
draw(upset_plot)

# Page 2: ggplot panels b, c, d
grid.newpage()
composite <- (p_b / (p_c | p_d)) +
  plot_annotation(
    tag_levels = list(c("b", "c", "d"))
  )
print(composite)
message("[caption] Supplementary Figure: External Database Comparison")

dev.off()
cat("Saved composite:", composite_pdf, "\n")

# ---- Save companion CSV ----
cat("=== Saving companion CSV ===\n")

companion <- atlas[, .(
  human_symbol,
  bulk_logFC,
  bulk_padj,
  is_deg,
  in_kegg_nafld = human_symbol %in% kegg_genes,
  in_wp_nafld = human_symbol %in% wp_genes,
  in_opentargets = human_symbol %in% ot_genes,
  in_any_external = human_symbol %in% any_external,
  atlas_unique_deg = is_deg & !(human_symbol %in% any_external),
  has_coloc
)]

# Add Open Targets score
companion <- merge(companion, ot_scores[, .(gene_symbol, ot_score)],
                   by.x = "human_symbol", by.y = "gene_symbol", all.x = TRUE)

fwrite(companion, file.path(OUT_DIR, "database_comparison_companion.csv"))
cat("Saved companion CSV:", file.path(OUT_DIR, "database_comparison_companion.csv"), "\n")

# Also save the recovery summary
fwrite(recovery_dt, file.path(OUT_DIR, "database_recovery_summary.csv"))
cat("Saved recovery summary:", file.path(OUT_DIR, "database_recovery_summary.csv"), "\n")

# Save statistics as text file
writeLines(stats_lines, file.path(OUT_DIR, "database_comparison_statistics.txt"))
cat("Saved statistics:", file.path(OUT_DIR, "database_comparison_statistics.txt"), "\n")

cat("\n=== DONE ===\n")
