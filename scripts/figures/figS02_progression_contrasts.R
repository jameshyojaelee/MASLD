#!/usr/bin/env Rscript
# figS02_progression_contrasts.R — Supplementary Figure 2
# "All progression contrast details" (5-panel composite)
# Regenerated with padj<0.05 thresholds (was padj<0.1)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(fgsea)
  library(msigdbr)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

prog_dir <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                      "results/progression")
out_dir  <- FIGS02_DIR

out_file <- file.path(out_dir, "figS02_progression_contrasts.pdf")

# ---------------------------------------------------------------------------
# (a) DEG count bar chart — binary + adjacent contrasts at padj<0.05
# ---------------------------------------------------------------------------
cat("=== Panel (a): DEG count bar chart ===\n")

binary_sum  <- fread(file.path(prog_dir, "progression_contrast_summary.csv"))
trans_sum   <- fread(file.path(prog_dir, "transition_summary.csv"))

# Build unified table for bar chart
bar_binary <- binary_sum[, .(contrast, n_deg = n_deg_05, type = "Binary")]
# Clean contrast names
bar_binary[, label := gsub("^C[0-9]+_", "", contrast)]

# Only fibrosis transitions for panel (a)
bar_trans <- trans_sum[stage_type == "fibrosis", .(contrast = transition,
                                             n_deg = n_deg_05,
                                             type = "Adjacent")]
bar_trans[, label := gsub("_to_", " \u2192 ", contrast)]

bar_dt <- rbind(bar_binary, bar_trans)
# Sort descending within type
bar_dt[, label := factor(label, levels = bar_dt[order(type, -n_deg), label])]

pa <- ggplot(bar_dt, aes(x = label, y = n_deg, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::comma(n_deg)), vjust = -0.3, size = GEOM_TEXT_6PT) +
  scale_fill_manual(values = c(Binary = masld_colors$deg,
                               Adjacent = "#F57F17"),
                    name = "Contrast type") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "DEGs (padj < 0.05)") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
message("[caption] Panel a: DEG counts per contrast")

# ---------------------------------------------------------------------------
# (b) Pairwise correlation heatmap — binary contrasts
# ---------------------------------------------------------------------------
cat("=== Panel (b): Pairwise correlation heatmap ===\n")

binary_files <- list.files(prog_dir, pattern = "^c[0-9]+.*_dream\\.csv$",
                           full.names = TRUE)
# Exclude ordinal contrasts (c7a, c7b, c7c, c17) — they are ordinal, not binary
binary_files <- binary_files[!grepl("c7[abc]|c17", basename(binary_files))]

# Read all, keep gene + logFC + a contrast id.
# C2 swap (2026-06-08): the binary c*_dream.csv files became limma-voom and dropped
# the `contrast_id` column — derive the contrast id from the file name instead.
read_binary <- function(f) {
  d <- fread(f, select = c("gene", "logFC"))
  cid <- sub("_dream\\.csv$", "", basename(f))
  d[, .(gene, logFC, cid = cid)]
}
all_binary <- rbindlist(lapply(binary_files, read_binary))

# Wide matrix
wide <- dcast(all_binary, gene ~ cid, value.var = "logFC")
mat  <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$gene

# Remove rows with any NA
mat <- mat[complete.cases(mat), ]

# Spearman correlation
cor_mat <- cor(mat, method = "spearman")

# Clean names for display
clean_name <- function(x) gsub("^c[0-9]+_", "", x)
rownames(cor_mat) <- clean_name(rownames(cor_mat))
colnames(cor_mat) <- clean_name(colnames(cor_mat))

# Melt for ggplot
cor_dt <- as.data.table(reshape2::melt(cor_mat, varnames = c("x", "y"),
                                       value.name = "rho"))

pb <- ggplot(cor_dt, aes(x = x, y = y, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = GEOM_TEXT_6PT) +
  scale_fill_gradient2(low = masld_colors$down, mid = "white",
                       high = masld_colors$up, midpoint = 0,
                       limits = c(-1, 1), name = "Spearman \u03C1") +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6))
message("[caption] Panel b: Pairwise logFC correlation")

# ---------------------------------------------------------------------------
# (c) Top pathway enrichment per fibrosis transition — run fgsea on the fly
# ---------------------------------------------------------------------------
cat("=== Panel (c): Pathway enrichment per transition ===\n")

trans_dream <- fread(file.path(prog_dir, "transition_fib_dream_results.csv"))
# Strip Ensembl version suffixes for gene set matching
trans_dream[, gene_nover := sub("\\..*", "", gene)]

# Get Hallmark gene sets
h_sets <- msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(h_sets$ensembl_gene, h_sets$gs_name)

run_fgsea_transition <- function(tr, dt) {
  sub <- dt[transition == tr]
  ranks <- setNames(sub$t, sub$gene_nover)
  ranks <- ranks[!is.na(ranks)]
  res <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500,
               nPermSimple = 10000)
  res$transition <- tr
  # Drop leadingEdge list-column for rbindlist compatibility
  res[, leadingEdge := NULL]
  res
}

fib_transitions <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
fgsea_res <- rbindlist(lapply(fib_transitions, run_fgsea_transition,
                              dt = trans_dream))

cat("  fgsea results:", nrow(fgsea_res), "rows across",
    uniqueN(fgsea_res$transition), "transitions\n")

# Top 5 per transition by padj
fgsea_res[, pathway_short := gsub("^HALLMARK_", "", pathway)]
fgsea_res[, pathway_short := gsub("_", " ", pathway_short)]
fgsea_res[, pathway_short := tolower(pathway_short)]
# Title case
fgsea_res[, pathway_short := tools::toTitleCase(pathway_short)]

top5 <- fgsea_res[order(padj), head(.SD, 5), by = transition]
top5[, transition_label := gsub("_to_", " -> ", transition)]

if (nrow(top5) > 0) {
  pc <- ggplot(top5, aes(x = NES, y = reorder(pathway_short, NES),
                         size = size, color = -log10(padj + 1e-50))) +
    geom_point() +
    scale_color_gradient(low = "grey70", high = masld_colors$up,
                         name = "-log10(padj)") +
    scale_size_continuous(range = c(1, 3.5), name = "Gene set size") +
    facet_wrap(~ transition_label, scales = "free_y", ncol = 2) +
    labs(x = "Normalized enrichment score", y = NULL) +
    theme_masld() +
    theme(axis.text.y = element_text(size = 6),
          strip.text = element_text(size = 6))
  message("[caption] Panel c: Top Hallmark pathways per fibrosis transition")
} else {
  pc <- placeholder("No fgsea results — check gene ID matching")
}

# ---------------------------------------------------------------------------
# (d) Volcano plot grid — 4 fibrosis transitions
# ---------------------------------------------------------------------------
cat("=== Panel (d): Volcano plots ===\n")

# Build gene -> symbol lookup. C2 swap (2026-06-08): binary contrast files no longer
# carry a `symbol` column — derive symbols from the gene→symbol map instead.
sym_lookup <- add_symbols(fread(binary_files[1], select = "gene"), "gene")[, .(gene, symbol)]

trans_dream[, padj := adj.P.Val]
trans_dream <- merge(trans_dream, sym_lookup, by = "gene", all.x = TRUE)
trans_dream[is.na(symbol), symbol := gene]

volc_dt <- trans_dream[transition %in% fib_transitions]
cat("  Volcano data:", nrow(volc_dt), "rows,",
    uniqueN(volc_dt$transition), "transitions\n")

volc_dt[, sig := fifelse(padj < 0.05, "DEG", "NS")]
volc_dt[, transition_label := gsub("_to_", " -> ", transition)]
volc_dt[, transition_label := factor(transition_label,
                                     levels = c("F0 -> F1", "F1 -> F2",
                                                "F2 -> F3", "F3 -> F4"))]
volc_dt[, neg_log10p := -log10(P.Value)]

# Cap extreme -log10(p) for display
volc_dt[neg_log10p > 50, neg_log10p := 50]

# Top 10 by |t| per transition for labeling
top_genes <- volc_dt[order(-abs(t)), head(.SD, 10), by = transition]

pd <- ggplot(volc_dt, aes(x = logFC, y = neg_log10p, color = sig)) +
  rasterize_layer(geom_point(size = 0.2, alpha = 0.4)) +
  geom_text_repel(data = top_genes,
                  aes(label = symbol),
                  size = GEOM_TEXT_6PT, max.overlaps = 15,
                  segment.size = 0.2, color = "black",
                  min.segment.length = 0) +
  scale_color_manual(values = c(DEG = masld_colors$up, NS = masld_colors$ns),
                     name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3,
             color = "grey50") +
  facet_wrap(~ transition_label, ncol = 2, scales = "free") +
  labs(x = "logFC", y = expression(-log[10](P))) +
  theme_masld() +
  theme(legend.position = "bottom")
message("[caption] Panel d: Adjacent fibrosis transition volcanos")

# ---------------------------------------------------------------------------
# (e) Threshold comparison — grouped bars
# ---------------------------------------------------------------------------
cat("=== Panel (e): Threshold comparison ===\n")

# Combine binary + transition summaries.
# Canonical effect-size threshold migrated 2026-04-27 from |LFC|>0.3 to |LFC|>0.5
# (LOO-CV stability; ~1.41x fold change). Sensitivity columns (lfc02/lfc03)
# remain in the producer table for the supp sensitivity panels.
thresh_binary <- binary_sum[, .(contrast = gsub("^C[0-9]+_", "", contrast),
                                `padj < 0.1`  = n_deg_01,
                                `padj < 0.05` = n_deg_05,
                                `padj < 0.05 + |LFC| > 0.5` = n_deg_05_lfc05)]
thresh_trans  <- trans_sum[stage_type == "fibrosis",
                           .(contrast = gsub("_to_", " \u2192 ", transition),
                             `padj < 0.1`  = n_deg_01,
                             `padj < 0.05` = n_deg_05,
                             `padj < 0.05 + |LFC| > 0.5` = n_deg_05_lfc05)]

thresh_all <- rbind(thresh_binary, thresh_trans)
thresh_m <- melt(thresh_all, id.vars = "contrast",
                 variable.name = "Threshold", value.name = "n_deg")
# Order contrasts by padj<0.05 count descending
order_levels <- thresh_all[order(-`padj < 0.05`), contrast]
thresh_m[, contrast := factor(contrast, levels = order_levels)]

pe <- ggplot(thresh_m, aes(x = contrast, y = n_deg, fill = Threshold)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c(`padj < 0.1` = "#90CAF9",
                               `padj < 0.05` = masld_colors$deg,
                               `padj < 0.05 + |LFC| > 0.5` = masld_colors$fibrosis)) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Number of DEGs") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom",
        legend.text = element_text(size = 6))
message("[caption] Panel e: Two-tier threshold comparison")

# ---------------------------------------------------------------------------
# Composite assembly
# ---------------------------------------------------------------------------
cat("=== Assembling composite figure ===\n")

top_row    <- pa + pb + plot_layout(widths = c(1.2, 1))
middle_row <- pc
bottom_row <- pd + pe + plot_layout(widths = c(1.2, 1))

composite <- top_row / middle_row / bottom_row +
  plot_layout(heights = c(1, 1.2, 1.3)) +
  plot_annotation(tag_levels = "a")

save_fig_tall(composite, out_file,
              width = fig_full_width, height = 12)

cat("Saved:", out_file, "\n")
cat("Done.\n")
