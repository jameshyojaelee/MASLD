#!/usr/bin/env Rscript
# figS_pathway_programs.R — Supplementary figure: Missing Pathway Programs
# Scores pyroptosis, necroptosis, UPR branches, mito dynamics, Hippo/YAP-TAZ,
# and senescence via ssGSEA, then tests association with fibrosis stage.
# Output: FIGS02_DIR/figS_pathway_programs.pdf + companion CSV

suppressPackageStartupMessages({
  library(edgeR)
  library(GSVA)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(grid)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

cat("=== figS_pathway_programs.R ===\n")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
cat("[1] Loading merged DGEList ...\n")
dge <- readRDS(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds"))
cat("  Samples:", ncol(dge), " Genes:", nrow(dge), "\n")

# logCPM matrix (rows = Ensembl IDs)
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# Load unified metadata (has fibrosis_stage)
meta <- read.csv(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
  stringsAsFactors = FALSE)
cat("  Metadata rows:", nrow(meta), "\n")

# Restrict to samples in the DGEList
meta <- meta[meta$sample_id %in% colnames(logcpm), ]

# Restrict to samples with fibrosis stage
meta <- meta[!is.na(meta$fibrosis_stage), ]
cat("  Samples with fibrosis_stage:", nrow(meta), "\n")

# Create F0-F4 labels
meta$fibrosis_stage <- factor(paste0("F", meta$fibrosis_stage),
                              levels = c("F0","F1","F2","F3","F4"))
meta <- meta[!is.na(meta$fibrosis_stage), ]

# Subset logCPM to these samples
logcpm <- logcpm[, meta$sample_id]
cat("  Final N =", ncol(logcpm), "\n")

# Load dream results for panel (c) — has Ensembl→symbol mapping
dream_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv")
dream <- fread(dream_file)
cat("  Dream results:", nrow(dream), "genes\n")

# Build Ensembl→symbol map from dream results
ens2sym <- setNames(dream$symbol, dream$gene)
ens2sym <- ens2sym[!is.na(ens2sym) & ens2sym != ""]

# Convert logCPM rownames to gene symbols (use symbol where available)
logcpm_symbols <- logcpm
ensembl_ids <- rownames(logcpm_symbols)
new_names <- ifelse(ensembl_ids %in% names(ens2sym), ens2sym[ensembl_ids], ensembl_ids)

# Handle duplicates: keep highest mean expression
mean_expr <- rowMeans(logcpm_symbols)
dup_syms <- new_names[duplicated(new_names)]
if (length(dup_syms) > 0) {
  keep <- rep(TRUE, length(new_names))
  for (s in unique(dup_syms)) {
    idx <- which(new_names == s)
    best <- idx[which.max(mean_expr[idx])]
    keep[setdiff(idx, best)] <- FALSE
  }
  logcpm_symbols <- logcpm_symbols[keep, ]
  new_names <- new_names[keep]
}
rownames(logcpm_symbols) <- new_names
cat("  Mapped to", sum(new_names != ensembl_ids[seq_along(new_names)]),
    "gene symbols\n")

# ---------------------------------------------------------------------------
# 2. Define gene sets
# ---------------------------------------------------------------------------
gene_sets <- list(
  Pyroptosis        = c("NLRP3","CASP1","CASP4","CASP5","GSDMD","GSDME",
                        "IL1B","IL18","PYCARD","NEK7","TXNIP","P2RX7"),
  Necroptosis       = c("RIPK1","RIPK3","MLKL","TNFRSF1A","TRADD","FADD",
                        "CASP8","CFLAR","ZBP1","HMGB1"),
  UPR_IRE1          = c("ERN1","XBP1","DNAJB9","EDEM1","HERPUD1","HSPA5","PDIA4"),
  UPR_PERK          = c("EIF2AK3","ATF4","DDIT3","PPP1R15A","TRIB3","ASNS","SLC7A11"),
  UPR_ATF6          = c("ATF6","MBTPS1","MBTPS2","CALR","PDIA6","HSP90B1"),
  Mito_Fission      = c("DNM1L","FIS1","MFF","MIEF2","MIEF1"),
  Mito_Fusion       = c("MFN1","MFN2","OPA1","OMA1"),
  Hippo_YAP_TAZ     = c("CYR61","CTGF","ANKRD1","AMOTL2","BIRC5","AXL"),
  Senescence        = c("CDKN1A","CDKN2A","TP53","RB1","SERPINE1","IL6",
                        "CXCL8","MMP1","MMP3","GLB1")
)

# Check gene availability against symbol-mapped matrix
avail <- rownames(logcpm_symbols)
for (gs in names(gene_sets)) {
  found <- gene_sets[[gs]][gene_sets[[gs]] %in% avail]
  missing <- setdiff(gene_sets[[gs]], avail)
  cat(sprintf("  %s: %d/%d found", gs, length(found), length(gene_sets[[gs]])))
  if (length(missing) > 0) cat(sprintf(" (missing: %s)", paste(missing, collapse=",")))
  cat("\n")
  gene_sets[[gs]] <- found
}

# Remove gene sets with <3 genes
gene_sets <- gene_sets[sapply(gene_sets, length) >= 3]
cat("  Retained", length(gene_sets), "gene sets with >=3 genes\n")

if (length(gene_sets) == 0) stop("No gene sets have >=3 genes in the expression matrix")

# ---------------------------------------------------------------------------
# 3. Run ssGSEA
# ---------------------------------------------------------------------------
cat("[2] Running ssGSEA ...\n")

ssgsea_param <- ssgseaParam(exprData = logcpm_symbols, geneSets = gene_sets)
scores <- gsva(ssgsea_param, verbose = TRUE)

cat("  ssGSEA done:", nrow(scores), "pathways x", ncol(scores), "samples\n")

# Compute fission/fusion ratio
if (all(c("Mito_Fission","Mito_Fusion") %in% rownames(scores))) {
  ratio <- scores["Mito_Fission",] - scores["Mito_Fusion",]
  scores <- rbind(scores, Fission_Fusion_Ratio = ratio)
}

# Long format for plotting
scores_df <- as.data.frame(t(scores)) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta[, c("sample_id","fibrosis_stage")], by = "sample_id") %>%
  filter(!is.na(fibrosis_stage)) %>%
  pivot_longer(cols = -c(sample_id, fibrosis_stage),
               names_to = "pathway", values_to = "score")

# ---------------------------------------------------------------------------
# 4. Statistical tests (Kruskal-Wallis per pathway)
# ---------------------------------------------------------------------------
cat("[3] Running Kruskal-Wallis tests ...\n")

kw_results <- scores_df %>%
  group_by(pathway) %>%
  summarise(
    kw_stat = kruskal.test(score ~ fibrosis_stage)$statistic,
    kw_pval = kruskal.test(score ~ fibrosis_stage)$p.value,
    .groups = "drop"
  ) %>%
  mutate(kw_padj = p.adjust(kw_pval, method = "BH"),
         sig_label = case_when(
           kw_padj < 0.001 ~ "***",
           kw_padj < 0.01  ~ "**",
           kw_padj < 0.05  ~ "*",
           TRUE ~ "ns"
         ))

print(kw_results)

# Stage-level means for heatmap
stage_means <- scores_df %>%
  group_by(pathway, fibrosis_stage) %>%
  summarise(mean_score = mean(score), se = sd(score)/sqrt(n()), .groups = "drop")

# ---------------------------------------------------------------------------
# 5. Panel (a): Pathway activity heatmap
# ---------------------------------------------------------------------------
cat("[4] Building panels ...\n")

# ggplot heatmap
hm_df <- stage_means %>%
  left_join(kw_results[, c("pathway","sig_label","kw_padj")], by = "pathway")

# Z-score within pathway
hm_df <- hm_df %>%
  group_by(pathway) %>%
  mutate(zscore = (mean_score - mean(mean_score)) / sd(mean_score)) %>%
  ungroup()

# Order pathways by F4 z-score
pathway_order <- hm_df %>%
  filter(fibrosis_stage == "F4") %>%
  arrange(zscore) %>%
  pull(pathway)

hm_df$pathway <- factor(hm_df$pathway, levels = pathway_order)

# Annotation: add sig stars at F4 position
sig_df <- hm_df %>%
  filter(fibrosis_stage == "F4") %>%
  select(pathway, sig_label) %>%
  distinct()

pa <- ggplot(hm_df, aes(x = fibrosis_stage, y = pathway, fill = zscore)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(data = sig_df,
            aes(x = "F4", y = pathway, label = sig_label),
            inherit.aes = FALSE, size = 2.5, hjust = -0.3) +
  scale_fill_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$up,
                       midpoint = 0, name = "Z-score") +
  labs(x = "Fibrosis stage", y = NULL, title = "Pathway activity across fibrosis") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6))

# ---------------------------------------------------------------------------
# 6. Panel (b): Key pathway trajectories
# ---------------------------------------------------------------------------
traj_pathways <- c("Pyroptosis","Necroptosis","UPR_IRE1","UPR_PERK","UPR_ATF6",
                   "Fission_Fusion_Ratio")
traj_pathways <- traj_pathways[traj_pathways %in% unique(scores_df$pathway)]

traj_df <- stage_means %>%
  filter(pathway %in% traj_pathways) %>%
  mutate(pathway = factor(pathway, levels = traj_pathways),
         stage_num = as.numeric(fibrosis_stage))

# Palette for trajectories
traj_cols <- c(
  Pyroptosis = "#C2185B",
  Necroptosis = "#880E4F",
  UPR_IRE1 = "#1565C0",
  UPR_PERK = "#42A5F5",
  UPR_ATF6 = "#64B5F6",
  Fission_Fusion_Ratio = "#7B1FA2"
)
traj_cols <- traj_cols[names(traj_cols) %in% traj_pathways]

pb <- ggplot(traj_df, aes(x = fibrosis_stage, y = mean_score,
                           color = pathway, group = pathway)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                width = 0.2, linewidth = 0.3) +
  scale_color_manual(values = traj_cols, name = "Pathway") +
  labs(x = "Fibrosis stage", y = "ssGSEA score (mean +/- SE)",
       title = "Pathway trajectories across fibrosis") +
  theme_masld() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"))

# ---------------------------------------------------------------------------
# 7. Panel (c): Individual gene DEG status (pyroptosis + necroptosis)
# ---------------------------------------------------------------------------
# Collect all pyroptosis + necroptosis genes (using original curated lists
# before filtering to available genes, so we show missing genes too)
pyro_genes <- c("NLRP3","CASP1","CASP4","CASP5","GSDMD","GSDME",
                "IL1B","IL18","PYCARD","NEK7","TXNIP","P2RX7")
necro_genes <- c("RIPK1","RIPK3","MLKL","TNFRSF1A","TRADD","FADD",
                 "CASP8","CFLAR","ZBP1","HMGB1")
deg_genes <- unique(c(pyro_genes, necro_genes))

# Match to dream results (using 'symbol' column)
dream_sub <- dream %>%
  filter(symbol %in% deg_genes) %>%
  select(symbol, logFC, padj) %>%
  mutate(
    sig = ifelse(!is.na(padj) & padj < 0.05, TRUE, FALSE),
    pathway_group = case_when(
      symbol %in% pyro_genes & symbol %in% necro_genes ~ "Both",
      symbol %in% pyro_genes ~ "Pyroptosis",
      symbol %in% necro_genes ~ "Necroptosis"
    )
  )

cat("  Panel (c): matched", nrow(dream_sub), "of", length(deg_genes),
    "pyroptosis+necroptosis genes in dream results\n")

# Order genes by logFC within pathway group
dream_sub <- dream_sub %>%
  arrange(pathway_group, logFC) %>%
  mutate(symbol = factor(symbol, levels = unique(symbol)))

if (nrow(dream_sub) > 0) {
  pc <- ggplot(dream_sub, aes(x = symbol, y = logFC, fill = logFC)) +
    geom_col(width = 0.7) +
    geom_point(data = dream_sub %>% filter(sig),
               aes(x = symbol, y = logFC),
               shape = 8, size = 1.2, color = "black", inherit.aes = FALSE) +
    scale_fill_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$up,
                         midpoint = 0, name = "logFC") +
    facet_wrap(~pathway_group, scales = "free_x", nrow = 1) +
    labs(x = NULL, y = "Dream logFC (MASLD vs Control)",
         title = "Individual gene evidence (DEG status)",
         caption = "* = padj < 0.05") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          strip.text = element_text(size = 6.5))
} else {
  pc <- placeholder("No pyroptosis/necroptosis genes\nfound in dream results")
}

# ---------------------------------------------------------------------------
# 8. Panel (d): Fission/fusion ratio boxplot
# ---------------------------------------------------------------------------
if ("Fission_Fusion_Ratio" %in% unique(scores_df$pathway)) {
  ratio_df <- scores_df %>%
    filter(pathway == "Fission_Fusion_Ratio")

  ratio_kw <- kw_results %>% filter(pathway == "Fission_Fusion_Ratio")
  pval_label <- sprintf("KW p = %.1e", ratio_kw$kw_pval)

  pd <- ggplot(ratio_df, aes(x = fibrosis_stage, y = score, fill = fibrosis_stage)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
    scale_fill_manual(values = fibrosis_stage_colors, guide = "none") +
    annotate("text", x = 3, y = max(ratio_df$score, na.rm = TRUE) * 1.05,
             label = pval_label, size = 2.5) +
    labs(x = "Fibrosis stage", y = "Fission - Fusion score",
         title = "Mitochondrial fission/fusion balance") +
    theme_masld()
} else {
  pd <- placeholder("Fission/fusion ratio\n(insufficient genes)")
}

# ---------------------------------------------------------------------------
# 9. Compose figure
# ---------------------------------------------------------------------------
cat("[5] Composing figure ...\n")

fig <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

out_dir <- FIGS02_DIR

out_pdf <- file.path(out_dir, "figS_pathway_programs.pdf")
save_fig_tall(fig, out_pdf, width = fig_full_width, height = 7)
cat("  Saved:", out_pdf, "\n")

# ---------------------------------------------------------------------------
# 10. Companion CSV
# ---------------------------------------------------------------------------
out_csv <- file.path(out_dir, "figS_pathway_programs_data.csv")
companion <- stage_means %>%
  left_join(kw_results, by = "pathway")
fwrite(companion, out_csv)
cat("  Saved:", out_csv, "\n")

out_gene_csv <- file.path(out_dir, "figS_pathway_programs_gene_deg.csv")
fwrite(dream_sub, out_gene_csv)
cat("  Saved:", out_gene_csv, "\n")

cat("=== DONE ===\n")
