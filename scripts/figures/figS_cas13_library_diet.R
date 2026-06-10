#!/usr/bin/env Rscript
# ===========================================================================
# Figure S_lib_2 (Mouse data landscape) + Figure S_lib_3 (Cross-diet replication)
# Cas13 library supplementary figures -- 6 CLEAN diet groups
# ===========================================================================
# KEY MESSAGE (S_lib_2): The mouse MASLD model landscape spans 6 mechanistically
# distinct diet groups covering >660 samples, with a shared core of high-confidence
# cross-diet DEGs that are concordant with human and genetically supported.
#
# KEY MESSAGE (S_lib_3): Cross-diet replication identifies a robust core of genes
# dysregulated across 5+ diet paradigms, with nutritional deficiency and metabolic
# excess families showing complementary but overlapping signatures.
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(UpSetR)
  library(grid)
  library(gridExtra)
  library(ggforce)
  library(scales)
})

# ---------------------------------------------------------------------------
# Source publication theme
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Output directory
OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load diet taxonomy
# ---------------------------------------------------------------------------
tax <- fread(file.path(BASE, "Cas13_Library_Design/data/diet_taxonomy.csv"))
diet_groups <- tax$diet_group
diet_labels <- tax$diet_label
diet_colors_vec <- setNames(tax$color_hex, tax$diet_group)
diet_order <- c("MCD", "CDAHFD", "Western", "HFD")  # NASH+Western merged 2026-05-29

# Short composition labels for x-axis (avoids text overflow)
diet_short_comp <- c(
  MCD     = "Met-choline\ndeficient",
  CDAHFD  = "Choline-def.\nHFD",
  Western = "HF + fructose\n± chol.",
  HFD     = "60% fat\nonly"
)

# Display names: Western = merged FPC/FFC/GAN + DIAMOND/WD family (2026-05-29).
diet_display <- c(
  MCD     = "MCD",
  CDAHFD  = "CDAHFD",
  Western = "Western",
  HFD     = "HFD"
)

# ---------------------------------------------------------------------------
# UP-DEG definition: ashr-shrunk effect size (2026-05-28). A fixed |logFC|
# cutoff is unfair across diets of different sample size (Western n=42 ...
# NASH n=249). ashr (M02c) shrinks logFC by its SE, so a uniform threshold on
# the shrunk effect is power-fair. lfsr<0.05 & shrunk_logFC>0.5 (matches the
# human dream Tier-1 magnitude).
# ---------------------------------------------------------------------------
LFSR_THR       <- 0.05
SHRUNK_LFC_THR <- 0.5

# ---------------------------------------------------------------------------
# Load per-diet DE results
# ---------------------------------------------------------------------------
de_dir <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
de_list <- lapply(diet_order, function(d) {
  f <- file.path(de_dir, paste0(d, "_de_results.csv"))
  dt <- fread(f)
  dt[, diet_group := d]
  dt[, gene_base := sub("\\..*", "", gene)]
  dt
})
names(de_list) <- diet_order

# ---------------------------------------------------------------------------
# Load mouse gene metadata (for biotype annotation)
# ---------------------------------------------------------------------------
mouse_meta <- fread(file.path(BASE,
  "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv"))
mouse_meta[, gene_base := mouse_ensembl_base]

# ---------------------------------------------------------------------------
# Load ortholog mapping (mouse Ensembl -> human symbol)
# ---------------------------------------------------------------------------
ortho <- fread(cmd = paste0("zcat ", file.path(BASE,
  "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz")))
ortho[, mouse_base := sub("\\..*", "", mouse_gene_ensembl)]
ortho[, human_base := sub("\\..*", "", human_gene_ensembl)]
ortho[, otype_rank := ifelse(ortholog_type == "ortholog_one2one", 1, 2)]
ortho <- ortho[order(mouse_base, otype_rank)][!duplicated(mouse_base)]

# ---------------------------------------------------------------------------
# Load multi-evidence atlas (for human DE concordance + COLOC)
# ---------------------------------------------------------------------------
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("human_symbol", "ensembl_id", "gene_biotype",
                           "bulk_logFC", "bulk_padj",
                           "coloc_susie_best_pp4", "coloc_abf_best_pp4"))
atlas[, human_base := sub("\\..*", "", ensembl_id)]
atlas[, is_human_deg := !is.na(bulk_padj) & bulk_padj < 0.05 & abs(bulk_logFC) > 0.5]
atlas[, has_coloc := (!is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0.5) |
                      (!is.na(coloc_abf_best_pp4) & coloc_abf_best_pp4 > 0.5)]

# ---------------------------------------------------------------------------
# Compute per-diet summary statistics
# ---------------------------------------------------------------------------
diet_summary <- rbindlist(lapply(diet_order, function(d) {
  dt <- de_list[[d]]
  n_total <- tax[diet_group == d, total_n]
  if (is.character(n_total)) n_total <- as.integer(gsub("~", "", n_total))
  n_up <- sum(dt$lfsr < LFSR_THR & dt$shrunk_logFC > SHRUNK_LFC_THR, na.rm = TRUE)
  n_down <- sum(dt$lfsr < LFSR_THR & dt$shrunk_logFC < -SHRUNK_LFC_THR, na.rm = TRUE)
  data.table(diet_group = d,
             diet_label = tax[diet_group == d, diet_label],
             composition = diet_short_comp[d],
             n_samples = n_total,
             n_degs = n_up + n_down,
             n_up = n_up,
             n_down = n_down)
}))
diet_summary[, diet_group := factor(diet_group, levels = diet_order)]

cat("=== Per-diet summary ===\n")
print(diet_summary)

# ===========================================================================
# Shared: identify upregulated DEGs per diet
# ===========================================================================
up_genes_per_diet <- lapply(diet_order, function(d) {
  de_list[[d]][lfsr < LFSR_THR & shrunk_logFC > SHRUNK_LFC_THR, gene_base]
})
names(up_genes_per_diet) <- diet_order

all_up_genes <- unique(unlist(up_genes_per_diet))
gene_diet_count <- data.table(gene_base = all_up_genes)
for (d in diet_order) {
  gene_diet_count[, (d) := gene_base %in% up_genes_per_diet[[d]]]
}
gene_diet_count[, n_diets := rowSums(.SD), .SDcols = diet_order]

lfc_dt <- rbindlist(lapply(diet_order, function(d) {
  de_list[[d]][, .(gene_base, logFC_abs = abs(shrunk_logFC),
                   logFC_val = shrunk_logFC, diet = d)]
}))
max_lfc <- lfc_dt[, .(max_lfc = max(logFC_abs, na.rm = TRUE)), by = gene_base]
gene_diet_count <- merge(gene_diet_count, max_lfc, by = "gene_base", all.x = TRUE)
gene_diet_count <- gene_diet_count[order(-n_diets, -max_lfc)]

# ===========================================================================
# FIGURE S_lib_2: MOUSE DATA LANDSCAPE
# ===========================================================================

# --- Panel A: Bar chart -- sample size + DEG count per diet group ---
# Use short label on x-axis
diet_summary[, x_label := paste0(diet_label, "\n", composition)]
diet_summary[, x_label := factor(x_label, levels = x_label)]

# Uniform 6 pt across all text for publication. geom_text size 2.1 mm ≈ 6 pt
# (size * ggplot2 .pt). Anchored to the diet-label size.
PANEL_FONT     <- 6
PANEL_GEOM_TXT <- 2.1

# Top: sample size
p_samples <- ggplot(diet_summary, aes(x = diet_group, y = n_samples, fill = diet_group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = diet_colors_vec, guide = "none") +
  geom_text(aes(label = n_samples), vjust = -0.3, size = PANEL_GEOM_TXT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_x_discrete(labels = diet_display) +
  labs(y = "N samples", x = NULL, title = "Samples per diet group") +
  theme_masld(base_size = 9) + theme_pub() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = PANEL_FONT),
        axis.title   = element_text(size = PANEL_FONT),
        plot.title   = element_text(size = PANEL_FONT, face = "bold"),
        plot.margin  = margin(2, 2, 0, 2))

# Bottom: UP DEGs, simple diet names on x-axis
p_degs <- ggplot(diet_summary, aes(x = diet_group, y = n_up, fill = diet_group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = diet_colors_vec, guide = "none") +
  geom_text(aes(label = comma(n_up)), vjust = -0.3, size = PANEL_GEOM_TXT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = comma) +
  scale_x_discrete(labels = diet_display) +
  labs(y = "UP DEGs\n(lfsr<0.05, shrunk LFC>0.5)", x = NULL) +
  theme_masld(base_size = 9) + theme_pub() +
  theme(axis.text.x = element_text(size = PANEL_FONT, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = PANEL_FONT),
        axis.title  = element_text(size = PANEL_FONT),
        plot.margin = margin(0, 2, 2, 2))

panel_A <- p_samples / p_degs + plot_layout(heights = c(1, 1.4))

# --- Panel B: Heatmap -- Per-gene x Per-diet significance ---
top50 <- gene_diet_count[1:min(50, nrow(gene_diet_count)), gene_base]

lfc_mat <- matrix(NA_real_, nrow = length(top50), ncol = length(diet_order),
                  dimnames = list(top50, diet_order))
lfsr_mat <- matrix(NA_real_, nrow = length(top50), ncol = length(diet_order),
                   dimnames = list(top50, diet_order))
for (d in diet_order) {
  dt <- de_list[[d]]
  idx <- match(top50, dt$gene_base)
  lfc_mat[, d] <- dt$shrunk_logFC[idx]   # ashr-shrunk effect size
  lfsr_mat[, d] <- dt$lfsr[idx]
}
lfc_mat_clipped <- pmin(pmax(lfc_mat, -3), 3)

# Map to gene symbols + annotation
symbol_map <- merge(
  data.table(gene_base = top50),
  mouse_meta[, .(gene_base, mouse_symbol_gtf, mouse_biotype)],
  by = "gene_base", all.x = TRUE)
symbol_map <- merge(symbol_map, ortho[, .(mouse_base, human_gene_symbol, human_base)],
                    by.x = "gene_base", by.y = "mouse_base", all.x = TRUE)
symbol_map <- merge(symbol_map, atlas[, .(human_base, is_human_deg, has_coloc)],
                    by.x = "human_base", by.y = "human_base", all.x = TRUE)
symbol_map <- symbol_map[match(top50, gene_base)]

row_labels <- ifelse(!is.na(symbol_map$mouse_symbol_gtf) & symbol_map$mouse_symbol_gtf != "",
                     symbol_map$mouse_symbol_gtf, symbol_map$gene_base)

biotype_simple <- ifelse(is.na(symbol_map$mouse_biotype), "other",
  ifelse(grepl("protein_coding", symbol_map$mouse_biotype), "PC",
    ifelse(grepl("lncRNA|lincRNA", symbol_map$mouse_biotype), "lncRNA", "other")))
human_conc <- ifelse(is.na(symbol_map$is_human_deg), "no ortholog",
  ifelse(symbol_map$is_human_deg, "concordant", "not DE"))
coloc_flag <- ifelse(is.na(symbol_map$has_coloc), "no data",
  ifelse(symbol_map$has_coloc, "PP4>0.5", "PP4<=0.5"))

biotype_cols <- c(PC = "#0D47A1", lncRNA = "#7B1FA2", other = "#9E9E9E")
human_cols <- c(concordant = "#00695C", "not DE" = "#F4A674", "no ortholog" = "#E0E0E0")
coloc_cols <- c("PP4>0.5" = "#C2185B", "PP4<=0.5" = "#BDBDBD", "no data" = "#E0E0E0")
col_labels <- diet_display[diet_order]

sig_marker <- function(j, i, x, y, width, height, fill) {
  if (!is.na(lfsr_mat[i, j]) && lfsr_mat[i, j] < LFSR_THR) {
    grid.circle(x = x, y = y, r = unit(0.8, "mm"),
                gp = gpar(fill = "black", col = NA))
  }
}

ht_col <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))

ha_row <- rowAnnotation(
  Biotype = biotype_simple,
  `Human DE` = human_conc,
  COLOC = coloc_flag,
  col = list(Biotype = biotype_cols, `Human DE` = human_cols, COLOC = coloc_cols),
  annotation_name_gp = gpar(fontsize = 5),
  annotation_legend_param = list(
    Biotype = list(title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 4)),
    `Human DE` = list(title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 4)),
    COLOC = list(title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 4))),
  simple_anno_size = unit(2.5, "mm"),
  gap = unit(0.3, "mm"))

ha_col <- HeatmapAnnotation(
  Diet = unname(diet_display[diet_order]),
  col = list(Diet = setNames(unname(diet_colors_vec[diet_order]),
                             unname(diet_display[diet_order]))),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Diet = list(title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 4))),
  simple_anno_size = unit(2.5, "mm"))

ht <- Heatmap(
  lfc_mat_clipped,
  name = "Shrunk logFC",
  col = ht_col,
  row_labels = row_labels,
  column_labels = col_labels,
  row_names_gp = gpar(fontsize = 4.5),
  column_names_gp = gpar(fontsize = 5.5),
  column_names_rot = 30,
  top_annotation = ha_col,
  right_annotation = ha_row,
  cell_fun = sig_marker,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  border = TRUE,
  heatmap_legend_param = list(
    title = "Shrunk logFC",
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 4),
    legend_height = unit(2, "cm")),
  width = unit(4, "cm"),
  height = unit(11, "cm"))

# --- Save S_lib_2: Panel A only ---
cat("Saving S_lib_2...\n")
pdf_file_2 <- file.path(OUT_DIR, "01_mouse_data_landscape.pdf")

panel_A_tagged <- panel_A +
  plot_annotation(tag_levels = list(c("A", ""))) &
  theme(plot.tag = element_text(size = 10, face = "bold"))
ggsave(pdf_file_2, panel_A_tagged,
       width = fig_half_width * 1.3, height = 3.5,
       device = if (capabilities("cairo")) cairo_pdf else pdf)
cat("  -> Saved:", pdf_file_2, "\n")


# ===========================================================================
# FIGURE S_lib_3: CROSS-DIET REPLICATION
# ===========================================================================

# --- Panel A: UpSet -- Cross-diet UP-DEG overlap ---
all_up <- unique(unlist(up_genes_per_diet))

upset_df <- data.frame(gene = all_up, stringsAsFactors = FALSE)
for (d in diet_order) {
  upset_df[[d]] <- as.integer(upset_df$gene %in% up_genes_per_diet[[d]])
}

gene_diet_n <- gene_diet_count[gene_base %in% all_up]
n_4plus <- sum(gene_diet_n$n_diets >= 4)
n_3plus <- sum(gene_diet_n$n_diets >= 3)
cat("UP DEGs in >=3 of 4 diets:", n_3plus, "\n")
cat("UP DEGs in all 4 diets:", n_4plus, "\n")

# --- Panel B: Proportional area diagram (2 diet families) ---
# LIDPAD (atherogenic) dropped 2026-05-28 — not disease-relevant for MASH library.
fam_nutdef <- unique(c(up_genes_per_diet[["MCD"]], up_genes_per_diet[["CDAHFD"]]))
fam_metexc <- unique(c(up_genes_per_diet[["Western"]], up_genes_per_diet[["HFD"]]))

n_nd <- length(fam_nutdef)
n_me <- length(fam_metexc)
n_nd_me <- length(intersect(fam_nutdef, fam_metexc))

cat("\n=== Diet family overlaps ===\n")
cat("Nutritional deficiency (MCD+CDAHFD):", n_nd, "\n")
cat("Metabolic excess (NASH+Western+HFD):", n_me, "\n")
cat("ND & ME (shared):", n_nd_me, "\n")

# Exclusive counts for annotation
nd_only <- n_nd - n_nd_me
me_only <- n_me - n_nd_me

max_n <- max(n_nd, n_me)
r_nd <- sqrt(n_nd / max_n) * 2.0
r_me <- sqrt(n_me / max_n) * 2.0

euler_df <- data.frame(
  x0 = c(-0.5, 0.9),
  y0 = c(0, 0),
  r = c(r_nd, r_me),
  family = c("Nutritional deficiency\n(MCD + CDAHFD)",
             "Metabolic excess\n(Western + HFD)"),
  stringsAsFactors = FALSE
)
euler_df$family <- factor(euler_df$family, levels = euler_df$family)

fam_colors <- c("Nutritional deficiency\n(MCD + CDAHFD)" = "#E64B35",
                "Metabolic excess\n(Western + HFD)" = "#00A087")

panel_B_euler <- ggplot() +
  geom_circle(data = euler_df, aes(x0 = x0, y0 = y0, r = r, fill = family),
              alpha = 0.30, linewidth = 0.4, color = "gray30") +
  scale_fill_manual(values = fam_colors, name = "Diet family") +
  # Exclusive counts
  annotate("text", x = -1.3, y = 0, label = comma(nd_only),
           size = 2.5, fontface = "bold") +
  annotate("text", x = 1.7, y = 0, label = comma(me_only),
           size = 2.5, fontface = "bold") +
  # Shared overlap
  annotate("text", x = 0.2, y = 0, label = comma(n_nd_me),
           size = 2.5, fontface = "bold", color = "#880E4F") +
  # Total counts per family (outside circles)
  annotate("text", x = -1.3, y = 1.6, label = paste0("N=", comma(n_nd)),
           size = 2.0, color = "#E64B35", fontface = "italic") +
  annotate("text", x = 1.7, y = 1.6, label = paste0("N=", comma(n_me)),
           size = 2.0, color = "#00A087", fontface = "italic") +
  coord_fixed(clip = "off") +
  labs(title = "Diet family UP-DEG overlap") +
  theme_void(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6, face = "bold"),
        plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
        plot.margin = margin(5, 5, 15, 5))

# --- Panel C: Cleveland dot chart -- Replication tier gene counts ---
all_up_with_biotype <- merge(
  gene_diet_count[gene_base %in% all_up],
  mouse_meta[, .(gene_base, mouse_biotype)],
  by = "gene_base", all.x = TRUE)
all_up_with_biotype[, biotype_simple := ifelse(
  grepl("protein_coding", mouse_biotype), "Protein-coding",
  ifelse(grepl("lncRNA|lincRNA", mouse_biotype), "lncRNA", "Other"))]

tier_counts <- all_up_with_biotype[biotype_simple %in% c("Protein-coding", "lncRNA"),
  .N, by = .(n_diets, biotype_simple)]

all_tiers <- CJ(n_diets = 1:4, biotype_simple = c("Protein-coding", "lncRNA"))
tier_counts <- merge(all_tiers, tier_counts, by = c("n_diets", "biotype_simple"), all.x = TRUE)
tier_counts[is.na(N), N := 0]
tier_counts[, n_diets_f := factor(n_diets)]

panel_C <- ggplot(tier_counts,
                  aes(x = N, y = n_diets_f, color = biotype_simple,
                      shape = biotype_simple)) +
  geom_segment(aes(x = 0, xend = N, y = n_diets_f, yend = n_diets_f),
               linewidth = 0.3, show.legend = FALSE) +
  geom_point(size = 2.5) +
  geom_text(aes(label = comma(N)), hjust = -0.3, size = 2.0, show.legend = FALSE) +
  scale_shape_manual(values = c("Protein-coding" = 16, lncRNA = 1), name = "Biotype") +
  scale_color_manual(values = c("Protein-coding" = "#0D47A1", lncRNA = "#7B1FA2"),
                     name = "Biotype") +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Gene count", y = "Replication tier\n(n diets UP, lfsr<0.05, shrunk LFC>0.5)",
       title = "Replication tier gene counts") +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

# --- Build UpSet-style plot in pure ggplot (avoids base/grid mixing issues) ---

# Compute intersection data for ggplot UpSet
# Binary membership matrix already exists as upset_df
bin_cols <- diet_order
bin_mat <- as.matrix(upset_df[, bin_cols])

# Create intersection signatures
upset_df$sig <- apply(bin_mat, 1, paste, collapse = "")
int_counts <- as.data.table(upset_df)[, .N, by = sig][order(-N)]

# Keep top 20 intersections
int_top <- int_counts[1:min(20, nrow(int_counts))]
int_top[, rank := .I]

# Decode signatures back to set membership
for (i in seq_along(diet_order)) {
  int_top[, (diet_order[i]) := as.integer(substr(sig, i, i))]
}

# Shared x-axis range so the top bars and bottom dot matrix line up exactly.
n_int <- nrow(int_top)
upset_xlim <- c(0.5, n_int + 0.5)
upset_xexp <- expansion(mult = c(0.01, 0.01))

# Intersection bar chart (top)
p_int_bars <- ggplot(int_top, aes(x = rank, y = N)) +
  geom_col(width = 0.6, fill = "gray30") +
  geom_text(aes(label = comma(N)), vjust = -0.3, size = 1.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = comma) +
  scale_x_continuous(limits = upset_xlim, expand = upset_xexp, breaks = NULL) +
  labs(y = "Intersection\nsize", x = NULL) +
  theme_masld() + theme_pub() +
  theme(plot.margin = margin(2, 2, 0, 2),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Dot matrix (bottom)
dot_data <- melt(int_top, id.vars = c("sig", "N", "rank"),
                 measure.vars = diet_order,
                 variable.name = "diet", value.name = "member")
dot_data[, diet := factor(diet, levels = rev(diet_order))]

# Connected segments for multi-set intersections
# Use numeric y for both dots and segments, then relabel
diet_levels_rev <- rev(diet_order)
dot_data[, y_num := match(as.character(diet), diet_levels_rev)]

seg_data <- dot_data[member == 1, .(y_min = min(y_num),
                                      y_max = max(y_num)),
                      by = rank]
seg_data <- seg_data[y_min != y_max]

p_dots <- ggplot(dot_data, aes(x = rank, y = y_num)) +
  geom_segment(data = seg_data, aes(x = rank, xend = rank,
                                      y = y_min, yend = y_max),
               linewidth = 0.5, color = "gray30", inherit.aes = FALSE) +
  geom_point(aes(color = ifelse(member == 1, "active", "inactive")),
             size = 1.8) +
  scale_color_manual(values = c(active = "gray20", inactive = "gray85"),
                     guide = "none") +
  scale_y_continuous(breaks = seq_along(diet_levels_rev),
                     labels = diet_levels_rev,
                     limits = c(0.3, length(diet_levels_rev) + 0.7),
                     expand = expansion(mult = 0)) +
  scale_x_continuous(limits = upset_xlim, expand = upset_xexp, breaks = NULL) +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(plot.margin = margin(0, 2, 2, 0),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank())

# Set size bar chart (left side) -- separate ggplot
set_sizes <- data.table(
  diet = diet_levels_rev,
  y_num = seq_along(diet_levels_rev),
  n = sapply(diet_levels_rev, function(d) length(up_genes_per_diet[[d]])))
set_colors_rev <- diet_colors_vec[diet_levels_rev]

# Horizontal bars: x = 0 to n, reversed so bars grow left
p_set_bars <- ggplot(set_sizes, aes(y = y_num, fill = diet)) +
  geom_rect(aes(xmin = 0, xmax = -n, ymin = y_num - 0.3, ymax = y_num + 0.3)) +
  geom_text(aes(x = -n, label = comma(n)), hjust = 1.1, size = 1.6) +
  scale_fill_manual(values = set_colors_rev, guide = "none") +
  scale_y_continuous(breaks = seq_along(diet_levels_rev),
                     labels = diet_display[diet_levels_rev],
                     limits = c(0.3, length(diet_levels_rev) + 0.7),
                     expand = expansion(mult = 0)) +
  scale_x_continuous(breaks = pretty(c(-max(set_sizes$n), 0), 3),
                     labels = function(x) comma(abs(x)),
                     expand = expansion(mult = c(0.30, 0))) +  # left headroom so the longest count label (e.g. 3,519) isn't clipped
  labs(x = "Set size", y = NULL) +
  theme_masld() + theme_pub() +
  theme(plot.margin = margin(0, 0, 2, 2),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = PUB_AXIS_TEXT, hjust = 1))

# Compose UpSet as 2x2 grid layout:
# top-left = empty, top-right = intersection bars
# bottom-left = set bars, bottom-right = dot matrix
p_empty <- ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

upset_gg <- (p_empty + p_int_bars + p_set_bars + p_dots) +
  plot_layout(ncol = 2, nrow = 2,
              widths = c(0.22, 0.78),
              heights = c(1.2, 1))

# --- Save S_lib_3 ---
cat("Saving S_lib_3...\n")
pdf_file_3 <- file.path(OUT_DIR, "02_cross_diet_replication.pdf")

# All panels are pure ggplot -- save directly with patchwork
panel_A_tagged <- upset_gg + plot_annotation(tag_levels = list("A")) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

panel_B_tagged <- panel_B_euler + labs(tag = "B") +
  theme(plot.tag = element_text(size = 8, face = "bold"))
panel_C_tagged <- panel_C + labs(tag = "C") +
  theme(plot.tag = element_text(size = 8, face = "bold"))
panel_BC <- (panel_B_tagged | panel_C_tagged) +
  plot_layout(widths = c(1.1, 0.9))

# Full figure: A on top, B+C on bottom
full_fig <- panel_A_tagged / panel_BC +
  plot_layout(heights = c(1.2, 1))

save_fig_tall(full_fig, pdf_file_3, width = fig_full_width, height = 9)

cat("  -> Saved:", pdf_file_3, "\n")
cat("\nDone. Output files:\n")
cat("  ", pdf_file_2, "\n")
cat("  ", pdf_file_3, "\n")
