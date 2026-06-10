#!/usr/bin/env Rscript
# =============================================================================
# figS_transition_tf_activity.R
# Per-Transition TF Activity — supports mid-stage (F1-F3) inflection narrative
#
# Approach: Use gene-level t-statistics from transition dream results,
# run decoupleR::run_wmean() with DoRothEA regulons (A+B+C) per transition
# to infer TF activity scores.
#
# Output:
#   FIGS02_DIR/figS_transition_tf_activity.pdf
#   RNA-seq/results/stratified_causal/transition_tf_activity.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(decoupleR)
  library(dorothea)
  library(reshape2)
})

# Paths
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

TRANS_FILE <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/transition_fib_dream_results.csv")
DREAM_FILE <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
OUT_FIG    <- file.path(FIGS02_DIR, "figS_transition_tf_activity.pdf")
OUT_CSV    <- file.path(BASE, "RNA-seq/results/stratified_causal/transition_tf_activity.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

cat("Loading transition dream results...\n")
trans <- fread(TRANS_FILE)

# Build Ensembl -> symbol mapping from main dream results
cat("Building gene ID mapping...\n")
dream_map <- fread(DREAM_FILE, select = c("gene", "symbol"))
dream_map <- dream_map[!is.na(symbol) & symbol != ""]
dream_map <- unique(dream_map)

# Merge symbol into transition results
trans <- merge(trans, dream_map, by = "gene", all.x = TRUE)
trans <- trans[!is.na(symbol) & symbol != ""]
# If duplicated symbols per transition, keep the one with highest absolute t
trans <- trans[order(-abs(t))][!duplicated(paste0(symbol, "_", transition))]

cat("Genes with symbols per transition:\n")
print(trans[, .N, by = transition])

# ------ DoRothEA regulons (A+B+C) ------
cat("Loading DoRothEA regulons...\n")
regulons <- dorothea_hs
regulons <- regulons[regulons$confidence %in% c("A", "B", "C"), ]
regulons <- regulons[, c("tf", "target", "mor")]
regulons <- as.data.frame(regulons)
cat(sprintf("  %d TF-target pairs, %d unique TFs\n", nrow(regulons), length(unique(regulons$tf))))

# ------ Run decoupleR per transition ------
transitions <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
tf_results_list <- list()

for (tr in transitions) {
  cat(sprintf("Running decoupleR for %s...\n", tr))
  sub <- trans[transition == tr]

  # Build a named t-statistic vector -> matrix (1 row = 1 "condition")
  tstat_vec <- setNames(sub$t, sub$symbol)

  # decoupleR expects a matrix: rows = genes, cols = samples/conditions
  mat <- matrix(tstat_vec, ncol = 1, dimnames = list(names(tstat_vec), tr))

  # Run weighted mean
  res <- run_wmean(mat = mat, net = regulons, .source = "tf", .target = "target",
                   .mor = "mor", times = 1000, minsize = 5)

  # Extract consensus/wmean results
  res_wmean <- res[res$statistic == "norm_wmean", ]
  res_wmean$transition <- tr
  tf_results_list[[tr]] <- as.data.table(res_wmean)
}

tf_all <- rbindlist(tf_results_list)
setnames(tf_all, "source", "tf")
setnames(tf_all, "condition", "condition_label")

cat(sprintf("Total TF-transition entries: %d\n", nrow(tf_all)))
cat(sprintf("Unique TFs: %d\n", length(unique(tf_all$tf))))

# Save full results
fwrite(tf_all[, .(tf, transition, score, p_value, statistic)], OUT_CSV)
cat(sprintf("Saved: %s\n", OUT_CSV))

# ------ Panel (a): TF activity heatmap — top 30 TFs ------
cat("Generating panel (a): heatmap...\n")

# Rank TFs by max |score| across transitions
tf_max <- tf_all[, .(max_abs = max(abs(score))), by = tf][order(-max_abs)]
top30 <- tf_max$tf[1:min(30, nrow(tf_max))]

# Build matrix
hm_data <- tf_all[tf %in% top30, .(tf, transition, score)]
hm_wide <- data.table::dcast(hm_data, tf ~ transition, value.var = "score")
hm_wide <- as.data.frame(hm_wide)
hm_mat <- as.matrix(hm_wide[, -1])
rownames(hm_mat) <- hm_wide$tf
# Order columns
hm_mat <- hm_mat[, transitions]

# Highlight key TFs
key_tfs <- c("RELA", "NFKB1", "NFKB2", "REL", "RELB", "SMAD3", "HNF4A", "PPARA", "STAT3", "ESR1")
key_present <- intersect(key_tfs, rownames(hm_mat))

# key_present used later for bold labels in ggplot heatmap
key_present <- intersect(key_tfs, rownames(hm_mat))

# ------ Panel (b): Key TF trajectory line plots ------
cat("Generating panel (b): trajectory line plots...\n")

plot_tfs <- c("RELA", "NFKB1", "SMAD3", "HNF4A", "PPARA", "STAT3", "ESR1", "NFKB2")
plot_tfs <- intersect(plot_tfs, unique(tf_all$tf))

line_data <- tf_all[tf %in% plot_tfs]
# Create ordered factor for transitions
line_data[, transition := factor(transition, levels = transitions)]
# Create a stage position for x-axis (midpoint of transition)
line_data[, stage_pos := as.numeric(transition)]
# Nice labels
line_data[, trans_label := factor(
  c("F0_to_F1" = "F0\u2192F1", "F1_to_F2" = "F1\u2192F2",
    "F2_to_F3" = "F2\u2192F3", "F3_to_F4" = "F3\u2192F4")[as.character(transition)],
  levels = c("F0\u2192F1", "F1\u2192F2", "F2\u2192F3", "F3\u2192F4")
)]

# Color TFs by functional group
tf_group <- c(
  RELA = "NF-\u03BAB", NFKB1 = "NF-\u03BAB", NFKB2 = "NF-\u03BAB",
  SMAD3 = "Fibrotic", HNF4A = "Metabolic", PPARA = "Metabolic",
  STAT3 = "Inflammatory", ESR1 = "Hormonal"
)
line_data[, tf_group := tf_group[tf]]

# Assign colors
group_colors <- c(
  "NF-\u03BAB" = "#C2185B",
  "Fibrotic" = "#7B1FA2",
  "Metabolic" = "#1565C0",
  "Inflammatory" = "#E91E63",
  "Hormonal" = "#AD1457"
)

pb <- ggplot(line_data, aes(x = trans_label, y = score, color = tf, group = tf)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  # Add vertical line at F2->F3 transition to highlight transition
  geom_vline(xintercept = 3, linetype = "dotted", color = "gray40", linewidth = 0.3) +
  annotate("text", x = 3.05, y = max(line_data$score) * 0.95,
           label = "mid-stage (F1-F3) inflection", hjust = 0, size = 2, color = "gray40", fontface = "italic") +
  scale_color_manual(values = setNames(
    c("#C2185B", "#E53935", "#D81B60",  # NF-kB family
      "#7B1FA2",                         # SMAD3
      "#1565C0", "#0D47A1",              # Metabolic
      "#E91E63",                         # STAT3
      "#AD1457"),                        # ESR1
    c("RELA", "NFKB1", "NFKB2", "SMAD3", "HNF4A", "PPARA", "STAT3", "ESR1")
  )) +
  labs(x = "Fibrosis transition", y = "TF activity score\n(normalized weighted mean)",
       color = "TF", title = "Key TF activity trajectories across fibrosis progression") +
  theme_masld(base_size = 7) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# ------ Panel (c): Volcano at F2->F3 ------
cat("Generating panel (c): F2->F3 volcano...\n")

f2f3 <- tf_all[transition == "F2_to_F3"]
f2f3[, neg_log10p := -log10(pmax(p_value, 1e-300))]

# Significance threshold
f2f3[, sig := ifelse(p_value < 0.05, ifelse(score > 0, "Activated", "Suppressed"), "NS")]

# Label top TFs
top_act <- f2f3[sig == "Activated"][order(-score)][1:min(8, sum(f2f3$sig == "Activated"))]
top_sup <- f2f3[sig == "Suppressed"][order(score)][1:min(8, sum(f2f3$sig == "Suppressed"))]
label_tfs <- unique(c(top_act$tf, top_sup$tf, intersect(key_tfs, f2f3$tf)))
f2f3[, label := ifelse(tf %in% label_tfs, tf, "")]

pc <- ggplot(f2f3, aes(x = score, y = neg_log10p, color = sig)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = c(
    "Activated" = masld_colors$up,
    "Suppressed" = masld_colors$down,
    "NS" = masld_colors$ns
  )) +
  ggrepel::geom_text_repel(
    data = f2f3[label != ""],
    aes(label = label),
    size = 2, max.overlaps = 20, segment.size = 0.2,
    min.segment.length = 0, fontface = "italic",
    color = "black"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "gray50") +
  labs(x = "TF activity score (F2\u2192F3)",
       y = expression(-log[10](p-value)),
       color = "Status",
       title = "TF activity at the F2\u2192F3 transition") +
  theme_masld(base_size = 7) +
  theme(legend.position = "right")

# ------ Composite figure ------
cat("Assembling composite figure...\n")

pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

# Use multi-page PDF: page 1 = heatmap, page 2 = ggplot panels (b+c)
# Then merge into single composite using a ggplot-based heatmap instead

# Convert heatmap to ggplot using geom_tile for clean patchwork integration
hm_gg_data <- as.data.table(melt(hm_mat))
setnames(hm_gg_data, c("tf", "transition", "score"))
hm_gg_data[, transition := factor(transition, levels = transitions)]
hm_gg_data[, trans_label := factor(
  c("F0_to_F1" = "F0\u2192F1", "F1_to_F2" = "F1\u2192F2",
    "F2_to_F3" = "F2\u2192F3", "F3_to_F4" = "F3\u2192F4")[as.character(transition)],
  levels = c("F0\u2192F1", "F1\u2192F2", "F2\u2192F3", "F3\u2192F4")
)]

# Cluster rows using same method as ComplexHeatmap
hc <- hclust(dist(hm_mat), method = "ward.D2")
tf_order <- hc$labels[hc$order]
hm_gg_data[, tf := factor(tf, levels = tf_order)]

# Bold key TFs
key_present <- intersect(key_tfs, levels(hm_gg_data$tf))
tf_face <- ifelse(levels(hm_gg_data$tf) %in% key_present, "bold.italic", "plain")

pa <- ggplot(hm_gg_data, aes(x = trans_label, y = tf, fill = score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                       midpoint = 0, name = "TF activity\n(norm. wmean)") +
  labs(x = NULL, y = NULL, title = "Per-transition TF activity (top 30)") +
  theme_masld(base_size = 7) +
  theme(axis.text.y = element_text(face = tf_face, size = 5.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3))

p_composite <- (pa | (pb / pc)) +
  plot_layout(widths = c(0.9, 1.1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 9, face = "bold", family = "Helvetica")))

save_fig(p_composite, OUT_FIG, width = fig_full_width, height = 5.5)
cat(sprintf("Saved figure: %s\n", OUT_FIG))

# ------ Summary stats ------
cat("\n=== Summary ===\n")
cat(sprintf("Total TFs tested: %d\n", length(unique(tf_all$tf))))
n_sig <- tf_all[p_value < 0.05, .N, by = transition]
cat("Significant TFs per transition (p < 0.05):\n")
print(n_sig)

# Key TFs at F2->F3
cat("\nKey TFs at F2->F3:\n")
f2f3_key <- tf_all[transition == "F2_to_F3" & tf %in% key_tfs, .(tf, score, p_value)][order(-score)]
print(f2f3_key)

cat("\nDone.\n")
