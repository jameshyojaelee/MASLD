#!/usr/bin/env Rscript
# figS04_coloc.R — Supplementary Figure 4: COLOC Results
# SuSiE-COLOC canonical (2026-07-05 MVP rebuild) across 50 GWAS strata
#   EUR 21 + AFR 10 + EAS 9 + AMR 7 + SAS 3 (27 MVP + 23 legacy studies)
# 18,975 eGenes tested.
# SuSiE PP.H4 > 0.5/0.8/0.9 = 736/569/416; ABF fallback 1,234/547.
# ARCHIVED 2026-04-09: Whitfield 36653562 removed (provenance unverified).
# Output: figures/supplementary/figS04_coloc/figS04_coloc.pdf (4-panel, 2x2)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

outdir <- FIGS04_DIR
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGS04_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
master <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
gene_level <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))

# Drop rows with empty gene names
master <- master[gene != "" & !is.na(gene)]
gene_level <- gene_level[gene != "" & !is.na(gene)]

cat("Master table:", nrow(master), "rows,", uniqueN(master$gwas_name), "GWAS\n")
cat("Gene-level table:", nrow(gene_level), "genes\n")

# ---------------------------------------------------------------------------
# GWAS phenotype category mapping
# ---------------------------------------------------------------------------
gwas_category <- data.table(
  gwas_name = c(
    "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
    "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
    "2021_34841290_NAFLD_EUR",
    "2023_36280732_NAFLD_deCode_EUR", "2023_36280732_NAFLD_Intermountain_EUR",
    "2023_36280732_NAFLD_UKBB_EUR",
    "FinnGen_NAFLD", "FinnGen_NASH",  # RESTORED 2026-04-09: verified FinnGen R12 provenance
    # ARCHIVED: Anstee_NAFLD (duplicate), 5x Whitfield 36653562 (unverified provenance), Pazoki_PDFF (duplicate)
    # DROPPED 2026-06: Ghouse_Cirrhosis, FinnGen_HCC, Ghouse_HCC (cirrhosis/HCC GWAS removed from 23-GWAS portfolio)
    "2021_34128465_PDFF_EUR", "2021_34957434_PDFF_EUR",
    "2022_36402844_PDFF_EUR"
  ),
  category = c(
    rep("Liver Enzymes", 3),
    rep("NAFLD/NASH", 8),
    rep("PDFF", 3)
  )
)

# Short labels for GWAS (biobank/method-based naming convention)
# 14 EUR GWAS after consolidation (7 archived, 3 cirrhosis/HCC dropped 2026-06)
gwas_labels <- data.table(
  gwas_name = gwas_category$gwas_name,
  gwas_short = c(
    "UKBB ALT", "UKBB AST", "UKBB GGT",
    "eMERGE NAFLD", "Biopsy NAFLD",
    "EHR NAFLD",
    "deCODE NAFLD", "Intermountain NAFLD",
    "UKBB NAFLD",
    "FinnGen NAFLD", "FinnGen NASH",
    "Abdominal MRI PDFF", "ML-derived PDFF",
    "Whole-body MRI PDFF"
  )
)

# Category colors (3 categories; Cirrhosis/HCC dropped with 23-GWAS portfolio)
category_colors <- c(
  "Liver Enzymes" = "#0D47A1",
  "NAFLD/NASH"    = "#C2185B",
  "PDFF"          = "#F57F17"
)

# ---------------------------------------------------------------------------
# Panel A: Per-GWAS hit counts (horizontal bar chart)
# ---------------------------------------------------------------------------
hits_per_gwas <- master[PP.H4.abf > 0.9, .N, by = gwas_name]
hits_per_gwas <- merge(hits_per_gwas, gwas_category, by = "gwas_name", all.x = TRUE)
hits_per_gwas <- merge(hits_per_gwas, gwas_labels, by = "gwas_name", all.x = TRUE)
hits_per_gwas[, gwas_short := factor(gwas_short, levels = gwas_short[order(N)])]

pA <- ggplot(hits_per_gwas, aes(x = N, y = gwas_short, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("n=", format(N, big.mark = ","))),
            hjust = -0.05, size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = category_colors, name = "Phenotype") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::comma) +
  labs(x = "Genes with PP.H4 > 0.9", y = NULL) +
  theme_masld(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        axis.text.y = element_text(size = 6))

# ---------------------------------------------------------------------------
# Panel B: Cross-GWAS replication histogram
# ---------------------------------------------------------------------------
repl_dist <- gene_level[, .(gene, n = coloc_n_gwas_h4_05)]
# Bin 6+ together
repl_dist[, n_bin := ifelse(n >= 6, "6+", as.character(n))]
repl_dist[, n_bin := factor(n_bin, levels = c("0", "1", "2", "3", "4", "5", "6+"))]

bin_counts <- repl_dist[, .N, by = n_bin]

# Top genes to annotate (4+ GWAS)
top_repl <- gene_level[coloc_n_gwas_h4_05 >= 4, .(gene, coloc_n_gwas_h4_05)]
top_repl <- top_repl[order(-coloc_n_gwas_h4_05)]

# Annotations for specific genes
annot_genes <- data.table(
  gene = c("C2orf16", "ATP13A1", "MTTP", "CWF19L1", "HKDC1", "EPHA2", "GPN1"),
  n_gwas = c(12, 8, 7, 7, 6, 6, 6)
)
annot_genes[, n_bin := ifelse(n_gwas >= 6, "6+", as.character(n_gwas))]

# Focus on genes with >=1 replication (drop the 0 bar which dominates)
bin_counts_nozero <- bin_counts[n_bin != "0"]

# Build annotation label for 6+ bar
annot_label <- paste(c("C2orf16 (12)", "ATP13A1 (8)", "MTTP (7)",
                        "CWF19L1 (7)", "HKDC1 (6)", "EPHA2 (6)"),
                      collapse = "\n")

message("[caption] Cross-GWAS replication: ",
        paste0(sum(gene_level$coloc_n_gwas_h4_05 >= 1, na.rm = TRUE),
               " genes colocalize in \u22651 GWAS; ",
               bin_counts[n_bin == "0", format(N, big.mark=",")],
               " with 0 omitted"))

pB <- ggplot(bin_counts_nozero, aes(x = n_bin, y = N)) +
  geom_col(fill = masld_colors$up, width = 0.7) +
  geom_text(aes(label = format(N, big.mark = ",")),
            vjust = -0.3, size = GEOM_TEXT_6PT, color = "black") +
  annotate("text", x = "5", y = bin_counts_nozero[n_bin == "1", N] * 0.85,
           label = annot_label,
           size = GEOM_TEXT_6PT, hjust = 0.5, color = "gray20", fontface = "italic",
           lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)),
                     labels = scales::comma) +
  labs(x = "Number of GWAS with PP.H4 > 0.5", y = "Number of genes") +
  theme_masld(base_size = 7)

# ---------------------------------------------------------------------------
# Panel C: Cross-phenotype heatmap
# ---------------------------------------------------------------------------
# Exclude MHC genes and LD-cluster secondary hits from top genes
gene_level_clean <- gene_level[is_mhc == FALSE | is.na(is_mhc)]
gene_level_clean <- gene_level_clean[is.na(ld_cluster_flag) |
                                      gene == ld_cluster_best_gene |
                                      is.na(ld_cluster_best_gene)]

# Select top genes: best by replication count, then by PP.H4
top_genes_repl <- gene_level_clean[coloc_n_gwas_h4_05 >= 2][order(-coloc_n_gwas_h4_05, -coloc_best_pp4)]
top_genes_pp4 <- gene_level_clean[coloc_best_pp4 >= 0.9][order(-coloc_best_pp4)]
top30_genes <- unique(c(top_genes_repl$gene, top_genes_pp4$gene))[1:30]

# Get per-gene per-category max PP.H4
master_cat <- merge(master, gwas_category, by = "gwas_name", all.x = TRUE)
heatmap_data <- master_cat[gene %in% top30_genes,
                           .(max_pp4 = max(PP.H4.abf, na.rm = TRUE)),
                           by = .(gene, category)]

# Complete grid
heatmap_grid <- CJ(gene = top30_genes,
                   category = names(category_colors))
heatmap_data <- merge(heatmap_grid, heatmap_data, by = c("gene", "category"), all.x = TRUE)
heatmap_data[is.na(max_pp4), max_pp4 := 0]

# Order genes by total replication
gene_order <- gene_level[gene %in% top30_genes][order(-coloc_n_gwas_h4_05, -coloc_best_pp4)]$gene
heatmap_data[, gene := factor(gene, levels = rev(gene_order))]
heatmap_data[, category := factor(category, levels = names(category_colors))]

# Known MASLD targets present in top 30
known_targets <- c("THRB", "SLC39A8", "DGAT2", "HSD17B13", "GCKR",
                    "PNPLA3", "TM6SF2", "NR1H4", "MBOAT7", "PPARA")
known_in_top30 <- intersect(known_targets, top30_genes)

# Mark known targets with a dot symbol in the label
gene_labels_c <- setNames(
  ifelse(rev(gene_order) %in% known_in_top30,
         paste0(rev(gene_order), " *"), rev(gene_order)),
  rev(gene_order)
)

message("[caption] Cross-phenotype colocalization (top 30 genes)")

pC <- ggplot(heatmap_data, aes(x = category, y = gene, fill = max_pp4)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colors = c("white", "#BBDEFB", "#1565C0", "#0D47A1"),
    values = scales::rescale(c(0, 0.3, 0.5, 1.0)),
    limits = c(0, 1),
    name = "Max PP.H4",
    breaks = c(0, 0.5, 0.8, 1.0)
  ) +
  scale_y_discrete(labels = gene_labels_c) +
  labs(x = NULL, y = NULL,
       caption = "* = known MASLD drug target; MHC-region and LD-cluster secondary genes excluded") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.25, "cm"),
        plot.caption = element_text(size = 6, color = "gray40", hjust = 0))

# ---------------------------------------------------------------------------
# Panel D: Known target validation (dot plot)
# ---------------------------------------------------------------------------
target_data <- gene_level[gene %in% known_targets,
                          .(gene, coloc_best_pp4, coloc_best_gwas,
                            coloc_n_gwas_h4_05)]
# Add category for best GWAS
target_data <- merge(target_data, gwas_category, by.x = "coloc_best_gwas",
                     by.y = "gwas_name", all.x = TRUE)

# Order by best PP.H4
target_data[, gene := factor(gene, levels = gene[order(coloc_best_pp4)])]

# Threshold line labels as data frames on the discrete scale
thresh_labels <- data.frame(
  gene = levels(target_data$gene)[nlevels(target_data$gene)],
  y = c(0.52, 0.82),
  label = c("PP.H4 = 0.9", "PP.H4 = 0.8"),
  color = c("gray50", "gray40")
)

pD <- ggplot(target_data, aes(x = gene, y = coloc_best_pp4)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40", linewidth = 0.3) +
  geom_text(data = thresh_labels, aes(x = gene, y = y, label = label),
            size = GEOM_TEXT_6PT, color = thresh_labels$color, hjust = 1.05, inherit.aes = FALSE) +
  geom_point(aes(color = category, size = coloc_n_gwas_h4_05)) +
  scale_color_manual(values = category_colors, name = "Best GWAS\ncategory") +
  scale_size_continuous(range = c(2, 5), name = "N GWAS\nPP.H4 > 0.9",
                        breaks = c(0, 1, 2, 3)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(x = "Known MASLD drug target", y = "Best PP.H4 (any GWAS)") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 6),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"))

# ---------------------------------------------------------------------------
# Assemble 2x2 figure
# ---------------------------------------------------------------------------
combined <- plot_grid(
  pA, pB,
  pC, pD,
  labels = c("a", "b", "c", "d"),
  label_size = 9,
  label_fontface = "plain",
  ncol = 2,
  rel_widths = c(1.1, 0.9),
  rel_heights = c(1, 1.1)
)

outfile <- file.path(outdir, "figS04_coloc.pdf")
save_fig_tall(combined, outfile, width = fig_full_width, height = 7.5)
cat("Saved:", outfile, "\n")

# ---------------------------------------------------------------------------
# Print summary statistics
# ---------------------------------------------------------------------------
cat("\n--- Summary ---\n")
cat("Total GWAS tested:", uniqueN(master$gwas_name), "\n")
cat("Genes with PP.H4 > 0.9 in any GWAS:", gene_level[coloc_n_gwas_h4_05 > 0, .N], "\n")
cat("Genes with PP.H4 > 0.8 in any GWAS:", gene_level[coloc_n_gwas_h4_08 > 0, .N], "\n")
cat("Genes replicating in 4+ GWAS:", gene_level[coloc_n_gwas_h4_05 >= 4, .N], "\n")
cat("Top replicated gene:", gene_level[which.max(coloc_n_gwas_h4_05), gene],
    "(", gene_level[which.max(coloc_n_gwas_h4_05), coloc_n_gwas_h4_05], "GWAS)\n")
cat("Known targets validated (PP.H4 > 0.9):",
    sum(target_data$coloc_best_pp4 > 0.5), "/", nrow(target_data), "\n")
cat("Known targets validated (PP.H4 > 0.8):",
    sum(target_data$coloc_best_pp4 > 0.8), "/", nrow(target_data), "\n")
