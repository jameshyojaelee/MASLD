#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: MASLD-Specific Network Analyses
#
# SBATCH: --partition=cpu, --mem=32G, --cpus-per-task=4, --time=48:00:00,
#         --job-name=figS_masld
#
# Addresses "hub genes are universal" critique and adds disease-specific
# insights to the multiplex gene network.
#
# Outputs:
#   A. Cross-layer posterior correlation matrix (7 x 7)
#   B. MASLD-specific hub genes (DEGs, COLOC, drug targets) + data CSV
#   C. Community x disease-category enrichment heatmap
#   D. mid-stage F1-F3 inflection gene network neighborhood
# ==========================================================================

t0 <- proc.time()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(ggraph)
  library(igraph)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
NET_DIR    <- NETWORK_DIR
COMM_DIR   <- file.path(NET_DIR, "communities")
OUT_DIR    <- FIGS_NET_DIR
DATA_DIR   <- file.path(OUT_DIR, "data")
dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# Guard: ID normalization must have run
id_norm_path <- file.path(NET_DIR, "id_normalization_summary.csv")
if (!file.exists(id_norm_path)) {
  stop("id_normalization_summary.csv not found -- ID normalization must run ",
       "before this figure script. Expected: ", id_norm_path)
}

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
ACTIVE_LAYERS <- c("ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna")
LAYER_COLS    <- paste0("p_", ACTIVE_LAYERS)

layer_labels <- c(
  ppi     = "PPI",
  coexpr  = "Co-expression",
  regulon = "Regulon",
  lr      = "Ligand-receptor",
  genetic = "Genetic",
  pathway = "Pathway",
  cerna   = "ceRNA"
)

COLOC_THR   <- 0.5            # coloc_susie_best_pp4 threshold
TOP_N_HUBS  <- 20
UNIVERSAL_HUBS <- c("ALB", "APOB", "APOE", "TTR", "FGB", "FGA", "SERPINA1",
                    "TP53", "AKT1", "EGFR", "MYC", "GAPDH", "ACTB")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
message("Loading nodes / composite edges / communities ...")
nodes <- fread(file.path(NET_DIR, "network_nodes.csv"))
setnames(nodes, "human_symbol", "gene")

comm <- fread(file.path(COMM_DIR, "community_assignments.csv"))
comm[, macro_id := as.character(macro_id)]

# Multi-evidence atlas for mid-stage (F1-F3) inflection gene identification
atlas_path <- file.path(ME, "multi_evidence_atlas.csv")
atlas <- if (file.exists(atlas_path)) fread(atlas_path) else NULL

# ===========================================================================
# A. Cross-layer posterior correlation matrix
# ===========================================================================
message("[A] computing cross-layer posterior correlations ...")

comp_path <- file.path(NET_DIR, "composite_edges.csv")
# Read only the layer columns (plus gene_a / gene_b for sanity)
comp_edges <- fread(comp_path, select = c("gene_a", "gene_b", LAYER_COLS))
# Keep only edges where at least 2 layers are non-zero, otherwise correlation
# is dominated by the mass of single-layer edges (which only carry info about
# their own layer).
n_layers_ok <- rowSums(comp_edges[, ..LAYER_COLS] > 0, na.rm = TRUE)
comp_multi <- comp_edges[n_layers_ok >= 2]
message(sprintf("  using %s edges with >=2 active layers (of %s total)",
                formatC(nrow(comp_multi), big.mark = ","),
                formatC(nrow(comp_edges), big.mark = ",")))

if (nrow(comp_multi) > 0) {
  # Pairwise spearman on edges where BOTH layers have values > 0
  L <- length(ACTIVE_LAYERS)
  cor_mat <- matrix(NA_real_, nrow = L, ncol = L,
                    dimnames = list(ACTIVE_LAYERS, ACTIVE_LAYERS))
  n_mat   <- matrix(NA_integer_, nrow = L, ncol = L,
                    dimnames = list(ACTIVE_LAYERS, ACTIVE_LAYERS))
  for (i in seq_len(L)) {
    for (j in seq_len(L)) {
      xi <- comp_multi[[LAYER_COLS[i]]]
      xj <- comp_multi[[LAYER_COLS[j]]]
      keep <- !is.na(xi) & !is.na(xj) & (xi > 0 | xj > 0)
      if (sum(keep) >= 50 && sd(xi[keep]) > 0 && sd(xj[keep]) > 0) {
        cor_mat[i, j] <- suppressWarnings(cor(xi[keep], xj[keep], method = "spearman"))
        n_mat[i, j]   <- sum(keep)
      }
    }
  }

  # Relabel for plotting
  rownames(cor_mat) <- colnames(cor_mat) <- layer_labels[ACTIVE_LAYERS]

  col_fun_corr <- colorRamp2(c(-0.5, 0, 0.5), c("#1565C0", "white", "#C2185B"))
  h_corr <- Heatmap(
    cor_mat,
    name = "Spearman\ncorrelation",
    col = col_fun_corr,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      v <- cor_mat[i, j]
      if (!is.na(v)) {
        grid.text(sprintf("%.2f", v), x, y,
                  gp = gpar(fontsize = 5,
                            col = ifelse(abs(v) > 0.35, "white", "black")))
      }
    },
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 6, fontface = "bold"),
      labels_gp = gpar(fontsize = 5),
      legend_height = unit(2.8, "cm")
    )
  )

  out_a <- file.path(OUT_DIR, "figS_masld_cross_layer_correlation.pdf")
  pdf(out_a, width = fig_col_width, height = 4.2)
  draw(h_corr,
       column_title = "Cross-layer posterior correlation",
       column_title_gp = gpar(fontsize = 8, fontface = "bold"))
  dev.off()
  message("  saved: ", out_a)

  fwrite(cbind(layer = rownames(cor_mat), as.data.frame(cor_mat)),
         file.path(DATA_DIR, "masld_cross_layer_correlation.csv"))
} else {
  message("  too few multi-layer edges to compute correlations")
}

rm(comp_multi); gc()

# ===========================================================================
# B. MASLD-specific hub genes
# ===========================================================================
# Composite degree restricted to each of 3 disease-relevant subsets.
# ===========================================================================
message("[B] computing MASLD-specific hubs ...")

# Free the full matrix before loading the lean edge table
rm(comp_edges); gc()

# Load only what we need: composite posterior + multiplicity
P_COMPOSITE_THR <- 0.5
hub_edges <- fread(comp_path,
                   select = c("gene_a", "gene_b", "p_composite", "k_multiplicity"))
hub_edges <- hub_edges[p_composite >= P_COMPOSITE_THR]
message(sprintf("  kept %s edges (p_composite >= %.2f)",
                formatC(nrow(hub_edges), big.mark = ","), P_COMPOSITE_THR))

compute_degree <- function(edges_dt, subset_genes) {
  eg <- edges_dt[gene_a %in% subset_genes & gene_b %in% subset_genes]
  if (nrow(eg) == 0) return(data.table(gene = character(), degree = integer()))
  dt <- rbind(
    eg[, .(gene = gene_a)],
    eg[, .(gene = gene_b)]
  )
  dt <- dt[, .(degree = .N), by = gene][order(-degree)]
  dt
}

# Subsets
deg_genes  <- nodes[!is.na(is_deg) & is_deg == TRUE, gene]
coloc_genes <- nodes[!is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > COLOC_THR, gene]
drug_genes  <- nodes[!is.na(dgidb_druggable) &
                     tolower(as.character(dgidb_druggable)) %in% c("yes", "true", "1"),
                     gene]

deg_hubs   <- compute_degree(hub_edges, deg_genes)
coloc_hubs <- compute_degree(hub_edges, coloc_genes)
drug_hubs  <- compute_degree(hub_edges, drug_genes)

# Also full-graph degree to compute universal-hub overlap
full_deg <- compute_degree(hub_edges,
                           unique(c(hub_edges$gene_a, hub_edges$gene_b)))
setorder(full_deg, -degree)
top_universal <- full_deg[1:min(TOP_N_HUBS, .N), gene]

top_deg    <- deg_hubs[1:min(TOP_N_HUBS, nrow(deg_hubs))]
top_coloc  <- coloc_hubs[1:min(TOP_N_HUBS, nrow(coloc_hubs))]
top_drug   <- drug_hubs[1:min(TOP_N_HUBS, nrow(drug_hubs))]

top_deg[,   subset := "DEG"]
top_coloc[, subset := "COLOC"]
top_drug[,  subset := "DrugTarget"]

# Annotate overlap with known universal hubs
mark_universal <- function(dt) {
  dt[, in_universal_top := gene %in% top_universal]
  dt[, in_known_universal := gene %in% UNIVERSAL_HUBS]
  dt
}
top_deg   <- mark_universal(top_deg)
top_coloc <- mark_universal(top_coloc)
top_drug  <- mark_universal(top_drug)

hub_table <- rbindlist(list(top_deg, top_coloc, top_drug))
hub_table[, rank := seq_len(.N), by = subset]

fwrite(hub_table, file.path(DATA_DIR, "masld_specific_hubs.csv"))
message("  data saved: ", file.path(DATA_DIR, "masld_specific_hubs.csv"))

# Plot: 3-panel bar chart, genes on y, degree on x, color by universal overlap
plot_hub_panel <- function(dt, title_str, bar_fill = "#C2185B") {
  dt <- copy(dt)
  dt[, gene := factor(gene, levels = rev(unique(gene)))]
  dt[, color_group := fifelse(in_known_universal, "Known universal hub",
                      fifelse(in_universal_top,   "Global hub (top 20)",
                                                  "Disease-specific"))]
  pal <- c("Disease-specific"     = bar_fill,
           "Global hub (top 20)"  = "#42A5F5",
           "Known universal hub"  = "#9E9E9E")
  ggplot(dt, aes(x = degree, y = gene, fill = color_group)) +
    geom_col(width = 0.75, color = "gray30", linewidth = 0.2) +
    scale_fill_manual(values = pal, drop = FALSE, name = NULL) +
    labs(title = title_str, x = "Composite degree (p>=0.5)", y = NULL) +
    theme_masld() +
    theme(axis.text.y = element_text(face = "italic", size = 6),
          axis.text.x = element_text(size = 6),
          plot.title  = element_text(face = "bold", size = 8),
          legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5))
}

p_b1 <- plot_hub_panel(top_deg,   sprintf("Top %d DEG hubs", TOP_N_HUBS),
                       bar_fill = "#C2185B")
p_b2 <- plot_hub_panel(top_coloc, sprintf("Top %d COLOC hubs (PP4>%.1f)",
                                          TOP_N_HUBS, COLOC_THR),
                       bar_fill = "#7B1FA2")
p_b3 <- plot_hub_panel(top_drug,  sprintf("Top %d druggable hubs", TOP_N_HUBS),
                       bar_fill = "#0D47A1")

p_b <- (p_b1 | p_b2 | p_b3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

out_b <- file.path(OUT_DIR, "figS_masld_specific_hubs.pdf")
save_fig(p_b, out_b, width = fig_full_width, height = 5.5)
message("  plot saved: ", out_b)

# ===========================================================================
# C. Community x disease-category enrichment heatmap
# ===========================================================================
message("[C] community x disease-category hypergeometric enrichment ...")

comm <- merge(comm, nodes[, .(gene, is_deg, coloc_susie_best_pp4,
                               dgidb_druggable, sex_class, ferroptosis_class,
                               zonation_class)], by = "gene", all.x = TRUE)

# Prepare boolean category vectors
categories <- list(
  DEGs          = comm$is_deg == TRUE,
  COLOC         = !is.na(comm$coloc_susie_best_pp4) & comm$coloc_susie_best_pp4 > COLOC_THR,
  DrugTargets   = !is.na(comm$dgidb_druggable) &
                   tolower(as.character(comm$dgidb_druggable)) %in% c("yes","true","1"),
  Sex_biased    = !is.na(comm$sex_class) & comm$sex_class %in% c("Male_biased", "Female_biased", "Divergent"),
  Ferroptosis   = !is.na(comm$ferroptosis_class) & nzchar(as.character(comm$ferroptosis_class)) &
                   comm$ferroptosis_class != "",
  Zonation      = !is.na(comm$zonation_class) & nzchar(as.character(comm$zonation_class)) &
                   comm$zonation_class != ""
)

comm_sizes <- comm[, .N, by = macro_id][order(-N)]
keep_macros <- comm_sizes[N >= 30, macro_id]

enrich_mat <- matrix(NA_real_, nrow = length(keep_macros), ncol = length(categories),
                     dimnames = list(keep_macros, names(categories)))
overlap_mat <- matrix(NA_integer_, nrow = length(keep_macros), ncol = length(categories),
                      dimnames = list(keep_macros, names(categories)))
or_mat <- matrix(NA_real_, nrow = length(keep_macros), ncol = length(categories),
                 dimnames = list(keep_macros, names(categories)))
N <- nrow(comm)

for (cat_i in seq_along(categories)) {
  cat_vec <- categories[[cat_i]]
  cat_vec[is.na(cat_vec)] <- FALSE
  K <- sum(cat_vec)
  if (K == 0) next
  for (m in keep_macros) {
    in_comm <- comm$macro_id == m
    n <- sum(in_comm)
    x <- sum(cat_vec & in_comm)
    if (x == 0 || n == 0) {
      enrich_mat[m, cat_i] <- 0
      overlap_mat[m, cat_i] <- x
      or_mat[m, cat_i] <- NA
      next
    }
    # Hypergeometric: P(X >= x)
    p <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
    enrich_mat[m, cat_i] <- -log10(pmax(p, 1e-300))
    # 2x2 OR
    a <- x; b <- n - x; c <- K - x; d <- N - K - (n - x)
    or <- tryCatch(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)),
                   error = function(e) NA)
    or_mat[m, cat_i] <- or
    overlap_mat[m, cat_i] <- x
  }
}

# BH correction across full matrix
padj_mat <- enrich_mat
p_vals <- 10 ^ -enrich_mat
p_vals_adj <- p.adjust(as.vector(p_vals), method = "BH")
padj_mat <- matrix(-log10(pmax(p_vals_adj, 1e-300)),
                   nrow = nrow(enrich_mat), ncol = ncol(enrich_mat),
                   dimnames = dimnames(enrich_mat))

# Order rows by size
row_order <- keep_macros
padj_mat <- padj_mat[row_order, , drop = FALSE]
or_mat   <- or_mat[row_order, , drop = FALSE]
overlap_mat <- overlap_mat[row_order, , drop = FALSE]

# Cap for visualization
plot_mat <- pmin(padj_mat, 50)

col_fun_enr <- colorRamp2(c(0, 1.3, 5, 15), c("white", "#FFEBEE", "#E91E63", "#880E4F"))

h_enr <- Heatmap(
  plot_mat,
  name = "-log10(padj)",
  col = col_fun_enr,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  row_title = "Macro community",
  column_title = "Disease category",
  row_title_gp = gpar(fontsize = 7, fontface = "bold"),
  column_title_gp = gpar(fontsize = 7, fontface = "bold"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- padj_mat[i, j]
    or <- or_mat[i, j]
    if (!is.na(v) && v >= 1.3) {
      label <- if (!is.na(or)) sprintf("%.1f\nOR=%.1f", v, or) else sprintf("%.1f", v)
      grid.text(label, x, y,
                gp = gpar(fontsize = 4.5,
                          col = ifelse(v > 8, "white", "black")))
    }
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5),
    legend_height = unit(2.8, "cm")
  )
)

out_c <- file.path(OUT_DIR, "figS_masld_community_enrichment.pdf")
pdf(out_c, width = fig_col_width, height = 5.0)
draw(h_enr,
     column_title = "Disease-category enrichment per community",
     column_title_gp = gpar(fontsize = 8, fontface = "bold"))
dev.off()
message("  saved: ", out_c)

fwrite(cbind(macro_id = rownames(padj_mat), as.data.frame(padj_mat)),
       file.path(DATA_DIR, "masld_community_enrichment_neglog10padj.csv"))
fwrite(cbind(macro_id = rownames(or_mat), as.data.frame(or_mat)),
       file.path(DATA_DIR, "masld_community_enrichment_OR.csv"))
fwrite(cbind(macro_id = rownames(overlap_mat), as.data.frame(overlap_mat)),
       file.path(DATA_DIR, "masld_community_enrichment_overlap.csv"))

# ===========================================================================
# D. mid-stage (F1-F3) inflection gene network neighborhood
# ===========================================================================
message("[D] mid-stage (F1-F3) inflection neighborhood ...")

# Identify mid-stage (F1-F3) inflection gene set.
# Preferred source: a precomputed list at results/convergence/ or similar.
# Fallback: inflammation-up + OxPhos-down DEG pattern from atlas.
F2_PATHWAY_UP   <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                     "HALLMARK_INFLAMMATORY_RESPONSE",
                     "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                     "HALLMARK_IL6_JAK_STAT3_SIGNALING")
F2_PATHWAY_DOWN <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                     "HALLMARK_FATTY_ACID_METABOLISM",
                     "HALLMARK_BILE_ACID_METABOLISM",
                     "HALLMARK_ADIPOGENESIS")

# Curated canonical MASLD mid-stage F1-F3 inflection hits (used when atlas-based selection
# returns too few genes).
F2_UP_FALLBACK   <- c("CCL2", "CCL20", "CXCL10", "IL6", "IL1B", "TNF",
                      "NFKB1", "RELA", "STAT1", "STAT3", "JUN", "FOS",
                      "COL1A1", "COL1A2", "COL3A1", "ACTA2", "TIMP1",
                      "TGFB1", "SERPINE1", "LGALS3", "SPP1", "MMP2")
F2_DOWN_FALLBACK <- c("PPARA", "CPT1A", "ACOX1", "HADHA", "HADHB",
                      "HMGCS2", "CYP2E1", "CYP7A1", "CYP1A2", "CYP3A4",
                      "ACADVL", "SLC27A5", "APOA5", "PCK1", "ALDOB")

# Prefer mid-stage (F1-F3) inflection gene list from previous analyses
f2_switch_candidates <- c(
  file.path(BASE, "RNA-seq/results/convergence/f2_switch_genes.csv"),
  file.path(BASE, "RNA-seq/results/convergence/switch_like_genes.csv"),
  file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/switch_like_genes.csv"),
  file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/f2_switch_genes.csv")
)
f2_file <- f2_switch_candidates[file.exists(f2_switch_candidates)][1]

f2_up   <- character(0)
f2_down <- character(0)
f2_source <- "none"

if (!is.na(f2_file)) {
  message("  using mid-stage (F1-F3) inflection list: ", f2_file)
  f2_dt <- fread(f2_file)
  sym_col <- intersect(c("symbol", "human_symbol", "gene", "hgnc_symbol"),
                       names(f2_dt))[1]
  dir_col <- intersect(c("direction", "switch_direction", "sign", "regulation"),
                      names(f2_dt))[1]
  if (!is.na(sym_col)) {
    genes_all <- as.character(f2_dt[[sym_col]])
    if (!is.na(dir_col)) {
      dir_vals <- as.character(f2_dt[[dir_col]])
      f2_up   <- genes_all[tolower(dir_vals) %in% c("up", "+", "positive", "1")]
      f2_down <- genes_all[tolower(dir_vals) %in% c("down", "-", "negative", "-1")]
    } else {
      f2_up <- genes_all
    }
    f2_source <- paste0("file:", basename(f2_file))
  }
}

# Atlas-based fallback
if ((length(f2_up) + length(f2_down)) < 20 && !is.null(atlas)) {
  message("  falling back to atlas-based inflammation-up + OxPhos-down pattern")
  sym_col <- intersect(c("human_symbol", "symbol"), names(atlas))[1]
  lfc_col <- intersect(c("bulk_logFC", "bulk_shrunk_logFC"), names(atlas))[1]
  padj_col <- intersect(c("bulk_padj", "bulk_lfsr"), names(atlas))[1]

  if (!is.na(sym_col) && !is.na(lfc_col) && !is.na(padj_col)) {
    atlas_slim <- atlas[!is.na(get(padj_col)) & get(padj_col) < 0.05 &
                        !is.na(get(lfc_col)) & abs(get(lfc_col)) > 0.5,
                        .(gene = get(sym_col),
                          lfc  = get(lfc_col),
                          padj = get(padj_col))]
    # Use pathway columns on atlas if present
    inflam_cols <- grep("HALLMARK_(TNFA|INFLAMM|IL6|INTERFERON_GAMMA)",
                         names(atlas), value = TRUE, ignore.case = TRUE)
    oxphos_cols <- grep("HALLMARK_(OXIDATIVE_PHOSPH|FATTY_ACID|BILE_ACID|ADIPOGEN)",
                         names(atlas), value = TRUE, ignore.case = TRUE)
    if (length(inflam_cols) > 0) {
      atlas[, f2_inflam := rowSums(as.matrix(.SD) > 0, na.rm = TRUE),
            .SDcols = inflam_cols]
      inflam_up_genes <- atlas[f2_inflam > 0 & !is.na(get(lfc_col)) & get(lfc_col) > 0,
                               get(sym_col)]
      f2_up <- unique(c(f2_up, inflam_up_genes))
    }
    if (length(oxphos_cols) > 0) {
      atlas[, f2_oxphos := rowSums(as.matrix(.SD) > 0, na.rm = TRUE),
            .SDcols = oxphos_cols]
      oxphos_dn_genes <- atlas[f2_oxphos > 0 & !is.na(get(lfc_col)) & get(lfc_col) < 0,
                               get(sym_col)]
      f2_down <- unique(c(f2_down, oxphos_dn_genes))
    }
    if (f2_source == "none") f2_source <- "atlas_pattern"
  }
}

# Curated fallback if atlas gave nothing
if ((length(f2_up) + length(f2_down)) < 10) {
  message("  using curated canonical mid-stage (F1-F3) inflection genes (fallback)")
  f2_up   <- unique(c(f2_up, F2_UP_FALLBACK))
  f2_down <- unique(c(f2_down, F2_DOWN_FALLBACK))
  if (f2_source == "none") f2_source <- "curated_fallback"
}

# Keep only genes present in network
f2_up   <- intersect(f2_up, nodes$gene)
f2_down <- intersect(f2_down, nodes$gene)

# Cap size for plotting readability
MAX_F2 <- 30
if (length(f2_up)   > MAX_F2) f2_up   <- f2_up[1:MAX_F2]
if (length(f2_down) > MAX_F2) f2_down <- f2_down[1:MAX_F2]
f2_genes <- unique(c(f2_up, f2_down))
message(sprintf("  F2 gene set (%s): %d up, %d down, %d unique",
                f2_source, length(f2_up), length(f2_down), length(f2_genes)))

# Extract edges among F2 genes + top bridge neighbors
f2_internal <- hub_edges[gene_a %in% f2_genes & gene_b %in% f2_genes &
                         p_composite >= 0.5]

# Add top external neighbors (p_composite >= 0.7) to surface community context
TOP_EXTERNAL <- 25
f2_bridge <- hub_edges[((gene_a %in% f2_genes) & !(gene_b %in% f2_genes) |
                        (gene_b %in% f2_genes) & !(gene_a %in% f2_genes)) &
                       p_composite >= 0.7]
setorder(f2_bridge, -p_composite)
if (nrow(f2_bridge) > TOP_EXTERNAL) f2_bridge <- f2_bridge[1:TOP_EXTERNAL]

f2_edges_all <- rbindlist(list(f2_internal, f2_bridge), fill = TRUE)
all_nodes_f2 <- unique(c(f2_edges_all$gene_a, f2_edges_all$gene_b))

if (length(all_nodes_f2) >= 2 && nrow(f2_edges_all) > 0) {
  node_meta <- data.frame(
    name = all_nodes_f2,
    group = ifelse(all_nodes_f2 %in% f2_up,   "F2 UP (inflam)",
            ifelse(all_nodes_f2 %in% f2_down, "F2 DOWN (OxPhos)", "Bridge")),
    stringsAsFactors = FALSE
  )
  # Degree within the subgraph
  g_f2 <- graph_from_data_frame(
    f2_edges_all[, .(from = gene_a, to = gene_b, p_composite, k_multiplicity)],
    vertices = node_meta,
    directed = FALSE
  )
  node_meta$degree <- degree(g_f2)
  V(g_f2)$degree <- node_meta$degree
  V(g_f2)$group  <- node_meta$group

  group_colors <- c(
    "F2 UP (inflam)"    = "#C2185B",
    "F2 DOWN (OxPhos)"  = "#1565C0",
    "Bridge"            = "#BDBDBD"
  )

  p_d <- ggraph(g_f2, layout = "fr") +
    geom_edge_link(aes(width = p_composite, alpha = p_composite),
                   color = "gray55") +
    geom_node_point(aes(size = degree, fill = group),
                    shape = 21, stroke = 0.3, color = "gray20") +
    geom_node_text(aes(label = name, fontface = ifelse(group == "Bridge", "plain", "bold")),
                   size = 1.8, repel = TRUE, max.overlaps = 40,
                   bg.color = "white", bg.r = 0.08) +
    scale_edge_width_continuous(range = c(0.15, 0.8), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
    scale_fill_manual(values = group_colors, name = "Role") +
    scale_size_continuous(range = c(1.5, 5), name = "Degree") +
    labs(title = "Mid-stage (F1-F3) metabolic-to-inflammatory inflection network neighborhood",
         subtitle = sprintf("%d nodes, %d edges  (source: %s; p_composite >= 0.5 internal, >= 0.7 bridge)",
                            vcount(g_f2), ecount(g_f2), f2_source)) +
    theme_masld() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), axis.title = element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 6, face = "bold"),
          plot.subtitle = element_text(size = 6, color = "gray40"),
          plot.title = element_text(face = "bold", size = 8))

  out_d <- file.path(OUT_DIR, "figS_midstage_inflection_neighborhood.pdf")
  save_fig(p_d, out_d, width = fig_full_width, height = 6.0)
  message("  saved: ", out_d)

  fwrite(data.table(gene = V(g_f2)$name, group = V(g_f2)$group,
                    degree = V(g_f2)$degree),
         file.path(DATA_DIR, "midstage_inflection_neighborhood_nodes.csv"))
  fwrite(f2_edges_all,
         file.path(DATA_DIR, "midstage_inflection_neighborhood_edges.csv"))
} else {
  message("  no F2 subgraph could be built (no edges found); skipping Panel D")
}

elapsed <- (proc.time() - t0)["elapsed"]
message(sprintf("figS_network_masld_specific.R completed in %.1f seconds", elapsed))
