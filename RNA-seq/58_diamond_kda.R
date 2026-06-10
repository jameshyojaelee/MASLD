#!/usr/bin/env Rscript
# =============================================================================
# Script 58: DIAMOnD Disease Module Expansion + Key Driver Analysis
# =============================================================================
# Gap #12: Network-based disease module detection and key driver identification
#
# DIAMOnD: Iteratively expands seed gene set (dream DEGs) using PPI connectivity
#   significance (hypergeometric test) to discover connector genes
# KDA: Identifies upstream key drivers in WGCNA modules via network neighborhood
#   enrichment in the STRING PPI
#
# Inputs:
#   - Bulk DEG results (canonical_deg_results.csv)
#   - WGCNA module assignments (wgcna_module_assignments.csv)
#   - STRING PPI (9606.protein.links.v12.0.txt.gz, score >= 700)
#   - SuSiE-COLOC gene-level results (gene_level_coloc.csv)
#   - Multi-evidence atlas (multi_evidence_atlas.csv)
#
# Outputs:
#   - figures/supplementary/figS08_subtyping_convergence/figS_diamond_kda.pdf (3 panels)
#   - RNA-seq/results/convergence/diamond_expanded_module.csv
#   - RNA-seq/results/convergence/key_driver_analysis.csv
#   - RNA-seq/results/convergence/diamond_kda_summary.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(igraph)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# --- Paths ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

dream_file   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
wgcna_file   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier/wgcna_module_assignments.csv")
string_links  <- file.path(BASE, "data/string_ppi/9606.protein.links.v12.0.txt.gz")
string_info   <- file.path(BASE, "data/string_ppi/9606.protein.info.v12.0.txt.gz")
coloc_file    <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
atlas_file    <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

out_fig    <- file.path(BASE, "figures/supplementary/figS08_subtyping_convergence/figS_diamond_kda.pdf")
out_dir    <- file.path(BASE, "RNA-seq/results/convergence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_fig), recursive = TRUE, showWarnings = FALSE)

cat("=== Script 58: DIAMOnD + KDA ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Load STRING PPI (score >= 700)
# =============================================================================
cat("Loading STRING PPI...\n")
ppi_raw <- fread(string_links)
setnames(ppi_raw, c("protein1", "protein2", "combined_score"))
ppi_raw <- ppi_raw[combined_score >= 700]
cat(sprintf("  High-confidence edges (score >= 700): %s\n", comma(nrow(ppi_raw))))

# Map ENSP -> gene symbol via STRING info file
info <- fread(string_info)
setnames(info, c("string_id", "preferred_name", "protein_size", "annotation"))
ensp2sym <- info[, .(string_id, symbol = preferred_name)]

ppi_raw[, sym1 := ensp2sym$symbol[match(protein1, ensp2sym$string_id)]]
ppi_raw[, sym2 := ensp2sym$symbol[match(protein2, ensp2sym$string_id)]]
ppi <- ppi_raw[!is.na(sym1) & !is.na(sym2) & sym1 != sym2]
cat(sprintf("  Mapped edges: %s\n", comma(nrow(ppi))))

# Build igraph
g <- graph_from_data_frame(ppi[, .(sym1, sym2)], directed = FALSE)
g <- simplify(g)
all_genes_ppi <- V(g)$name
cat(sprintf("  PPI network: %s nodes, %s edges\n",
            comma(vcount(g)), comma(ecount(g))))

# =============================================================================
# 2. Load Dream DEGs + COLOC genes
# =============================================================================
cat("\nLoading seed gene sets...\n")

# Dream DEGs
dream <- fread(dream_file)
# Need symbol mapping: dream uses ensembl IDs, atlas has the mapping
atlas <- fread(atlas_file, select = c("human_symbol", "ensembl_id",
                                       "bulk_logFC", "bulk_padj"))

# Primary DEGs
deg_all <- atlas[!is.na(bulk_padj) & bulk_padj < 0.05 & abs(bulk_logFC) > 0.3]
deg_symbols <- deg_all$human_symbol
deg_in_ppi <- intersect(deg_symbols, all_genes_ppi)
cat(sprintf("  DEGs (padj<0.05, |LFC|>0.3): %d total, %d in PPI\n",
            length(deg_symbols), length(deg_in_ppi)))

# Top 500 by |logFC| for focused seed set
top500 <- atlas[!is.na(bulk_padj) & bulk_padj < 0.05][
  order(-abs(bulk_logFC))][1:min(500, .N)]$human_symbol
top500_in_ppi <- intersect(top500, all_genes_ppi)
cat(sprintf("  Top-500 DEGs by |LFC|: %d in PPI\n", length(top500_in_ppi)))

# COLOC genes (PP.H4 > 0.5)
coloc <- fread(coloc_file)
coloc_genes <- coloc[coloc_best_pp4 > 0.5 & gene != ""]$gene
coloc_in_ppi <- intersect(coloc_genes, all_genes_ppi)
cat(sprintf("  COLOC genes (PP4>0.5): %d total, %d in PPI\n",
            length(coloc_genes), length(coloc_in_ppi)))

# =============================================================================
# 3. DIAMOnD Algorithm
# =============================================================================
# DIAMOnD: Disease Module Detection (Ghiassian et al., PLoS Comp Bio 2015)
# Iteratively adds the gene with the best connectivity significance to the
# current seed module.
cat("\nRunning DIAMOnD algorithm (top-500 seed)...\n")

# Precompute sparse adjacency matrix (shared between DIAMOnD and KDA)
cat("  Computing sparse adjacency matrix...\n")
A <- as_adjacency_matrix(g, sparse = TRUE)  # dgCMatrix
node_idx <- setNames(seq_along(all_genes_ppi), all_genes_ppi)
degree_vec <- Matrix::rowSums(A)
cat(sprintf("  Adjacency matrix: %d x %d\n", nrow(A), ncol(A)))

diamond <- function(seeds, max_added = 200) {
  N <- length(all_genes_ppi)
  seed_set <- intersect(seeds, all_genes_ppi)

  # Seed indicator vector
  seed_idx <- node_idx[seed_set]
  module_vec <- rep(0, N)
  module_vec[seed_idx] <- 1

  in_module <- rep(FALSE, N)
  in_module[seed_idx] <- TRUE

  results_list <- vector("list", max_added)

  for (i in seq_len(max_added)) {
    s0 <- sum(in_module)

    # k_vals = number of module neighbors for each node (matrix-vector product)
    k_vals <- as.numeric(A %*% module_vec)

    # Mask out genes already in module
    k_vals[in_module] <- NA

    deg_vals <- degree_vec
    deg_vals[in_module] <- NA

    # Hypergeometric p-value for all candidates at once
    p_vals <- phyper(k_vals - 1L, s0, N - s0, deg_vals, lower.tail = FALSE)
    p_vals[in_module] <- NA

    idx_best <- which.min(p_vals)
    if (length(idx_best) == 0) break

    best_gene <- all_genes_ppi[idx_best]

    in_module[idx_best] <- TRUE
    module_vec[idx_best] <- 1

    results_list[[i]] <- data.table(
      gene = best_gene,
      iteration = i,
      p_value = p_vals[idx_best],
      connectivity_to_module = as.integer(k_vals[idx_best]),
      degree = as.integer(degree_vec[idx_best]),
      module_size = sum(in_module)
    )

    if (i %% 50 == 0) {
      cat(sprintf("  Iteration %d: added %s (p=%.2e, k=%d)\n",
                  i, best_gene, p_vals[idx_best], k_vals[idx_best]))
    }
  }

  rbindlist(results_list[!vapply(results_list, is.null, logical(1))])
}

diamond_results <- diamond(top500_in_ppi, max_added = 200)
cat(sprintf("  DIAMOnD complete: %d connector genes added\n", nrow(diamond_results)))

# Annotate connector genes with DEG/COLOC status
diamond_results[, is_deg := gene %in% deg_symbols]
diamond_results[, is_coloc := gene %in% coloc_genes]
diamond_results[, category := fifelse(
  is_deg & is_coloc, "DEG+COLOC",
  fifelse(is_deg, "DEG_only",
  fifelse(is_coloc, "COLOC_only", "Connector"))
)]

# Save
fwrite(diamond_results, file.path(out_dir, "diamond_expanded_module.csv"))
cat(sprintf("  Connector category breakdown:\n"))
print(diamond_results[, .N, by = category][order(-N)])

# =============================================================================
# 4. Key Driver Analysis (KDA)
# =============================================================================
cat("\nRunning Key Driver Analysis...\n")

# Load WGCNA module assignments
wgcna <- fread(wgcna_file)
# Strip version from ensembl IDs to match atlas
wgcna[, ensembl_clean := sub("\\..*", "", gene)]

# Map to symbol via atlas
atlas_map <- atlas[, .(human_symbol, ensembl_clean = sub("\\..*", "", ensembl_id))]
wgcna <- merge(wgcna, atlas_map, by = "ensembl_clean", all.x = TRUE)
wgcna <- wgcna[!is.na(human_symbol) & human_symbol != ""]

# Remove grey module (unassigned)
wgcna_assigned <- wgcna[module_color != "grey"]
modules <- unique(wgcna_assigned$module_color)
cat(sprintf("  WGCNA modules (excluding grey): %d\n", length(modules)))
cat(sprintf("  Genes with module assignments in PPI: %d\n",
            sum(wgcna_assigned$human_symbol %in% all_genes_ppi)))

# KDA: For each gene in PPI, test enrichment of its neighbors in each module
# Use hypergeometric test (one-sided)

# Reuse sparse adjacency matrix A, node_idx, degree_vec from DIAMOnD
cat("  Reusing sparse adjacency matrix for KDA...\n")
N_kda <- length(all_genes_ppi)
# Candidates: genes with >= 3 neighbors
kda_cand_mask <- degree_vec >= 3
cat(sprintf("  KDA candidates (degree >= 3): %d\n", sum(kda_cand_mask)))

kda_run <- function(module_genes, all_ppi_genes, module_name) {
  module_in_ppi <- intersect(module_genes, all_ppi_genes)
  if (length(module_in_ppi) < 5) return(NULL)

  M <- length(module_in_ppi)
  N <- N_kda

  # Module indicator vector
  mod_vec <- rep(0, nrow(A))
  mod_vec[node_idx[module_in_ppi]] <- 1

  # k = number of module neighbors for each gene
  k_all <- as.numeric(A %*% mod_vec)

  # Filter: candidates with degree >= 3 and k >= 2
  keep <- kda_cand_mask & k_all >= 2
  if (sum(keep) == 0) return(NULL)

  idx_keep <- which(keep)
  k_f <- as.integer(k_all[idx_keep])
  n_f <- as.integer(degree_vec[idx_keep])
  cands <- all_genes_ppi[idx_keep]

  p_vals <- phyper(k_f - 1L, M, N - M, n_f, lower.tail = FALSE)

  data.table(
    gene = cands,
    module = module_name,
    module_neighbors = k_f,
    total_neighbors = n_f,
    module_size_in_ppi = M,
    fold_enrichment = (k_f / n_f) / (M / N),
    p_value = p_vals
  )
}

# Run KDA for each WGCNA module
kda_all <- rbindlist(lapply(modules, function(mod) {
  mod_genes <- wgcna_assigned[module_color == mod]$human_symbol
  cat(sprintf("  Module %s: %d genes (%d in PPI)\n",
              mod, length(mod_genes), sum(mod_genes %in% all_genes_ppi)))
  kda_run(mod_genes, all_genes_ppi, mod)
}))

if (nrow(kda_all) > 0) {
  # BH correction within each module
  kda_all[, padj := p.adjust(p_value, method = "BH"), by = module]
  kda_all <- kda_all[order(p_value)]

  # Annotate with DEG / COLOC status
  kda_all[, is_deg := gene %in% deg_symbols]
  kda_all[, is_coloc := gene %in% coloc_genes]

  # Save
  fwrite(kda_all, file.path(out_dir, "key_driver_analysis.csv"))

  # Summary: top drivers per module
  kda_top <- kda_all[padj < 0.05, .SD[which.min(p_value)], by = module]
  cat(sprintf("\n  Significant key drivers (padj < 0.05): %d genes across %d modules\n",
              kda_all[padj < 0.05, uniqueN(gene)],
              kda_all[padj < 0.05, uniqueN(module)]))
  cat("  Top driver per module:\n")
  print(kda_top[, .(module, gene, fold_enrichment = round(fold_enrichment, 2),
                     padj = signif(padj, 3), module_neighbors, is_deg, is_coloc)])
} else {
  cat("  WARNING: No KDA results generated\n")
}

# =============================================================================
# 5. Summary table
# =============================================================================
# Combine DIAMOnD + KDA highlights
summary_dt <- data.table(
  metric = c(
    "DIAMOnD_seed_genes_in_PPI",
    "DIAMOnD_connector_genes_added",
    "DIAMOnD_connectors_also_DEG",
    "DIAMOnD_connectors_also_COLOC",
    "DIAMOnD_connectors_novel",
    "DIAMOnD_sig_connectors_p005",
    "KDA_total_significant_drivers",
    "KDA_modules_with_drivers",
    "KDA_drivers_also_DEG",
    "KDA_drivers_also_COLOC"
  ),
  value = c(
    length(top500_in_ppi),
    nrow(diamond_results),
    sum(diamond_results$is_deg),
    sum(diamond_results$is_coloc),
    sum(diamond_results$category == "Connector"),
    sum(diamond_results$p_value < 0.05),
    if (nrow(kda_all) > 0) kda_all[padj < 0.05, uniqueN(gene)] else 0,
    if (nrow(kda_all) > 0) kda_all[padj < 0.05, uniqueN(module)] else 0,
    if (nrow(kda_all) > 0) kda_all[padj < 0.05 & is_deg == TRUE, uniqueN(gene)] else 0,
    if (nrow(kda_all) > 0) kda_all[padj < 0.05 & is_coloc == TRUE, uniqueN(gene)] else 0
  )
)
fwrite(summary_dt, file.path(out_dir, "diamond_kda_summary.csv"))
cat("\nSummary:\n")
print(summary_dt)

# =============================================================================
# 6. Figures (3 panels)
# =============================================================================
cat("\nGenerating figures...\n")

# --- Panel A: DIAMOnD expansion curve ---
pa <- ggplot(diamond_results, aes(x = iteration, y = module_size)) +
  geom_line(color = masld_colors$up, linewidth = 0.8) +
  geom_point(aes(color = -log10(p_value)), size = 1.2) +
  scale_color_gradient(low = "grey70", high = masld_colors$up,
                       name = expression(-log[10](p))) +
  geom_hline(yintercept = length(top500_in_ppi), linetype = "dashed",
             color = masld_colors$down, linewidth = 0.4) +
  annotate("text", x = 10, y = length(top500_in_ppi) + 15,
           label = paste0("Seed: ", length(top500_in_ppi), " DEGs"),
           color = masld_colors$down, hjust = 0, size = 3) +
  labs(x = "DIAMOnD Iteration",
       y = "Cumulative Module Size",
       title = "DIAMOnD Disease Module Expansion") +
  theme_masld()

# Inset: p-value trajectory
pa_inset <- ggplot(diamond_results, aes(x = iteration, y = -log10(p_value))) +
  geom_line(color = masld_colors$up, linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(x = "Iteration", y = expression(-log[10](p))) +
  theme_masld(base_size = 7) +
  theme(plot.background = element_rect(fill = "white", color = "grey80"))

# Wrap as single grob so inset does not consume a tag slot
pa_with_inset <- pa + inset_element(pa_inset, left = 0.55, bottom = 0.55,
                                     right = 0.98, top = 0.98)
pa_final <- wrap_elements(full = pa_with_inset)

# --- Panel B: Top 5 key drivers per module (diverse representation) ---
if (nrow(kda_all) > 0 && kda_all[padj < 0.05, .N] > 0) {
  # Select top 5 per module to ensure diversity across all WGCNA modules
  kda_top_per_mod <- kda_all[padj < 0.05][order(p_value), .SD[1:min(5, .N)], by = module]
  kda_top_per_mod[, gene := factor(gene, levels = rev(gene))]

  # Module-specific colors for visual grouping
  mod_pal <- c("blue" = "#1565C0", "brown" = "#8D6E63",
               "turquoise" = "#00897B", "yellow" = "#FDD835",
               "green" = "#43A047", "red" = "#E53935",
               "black" = "#424242", "pink" = "#EC407A",
               "magenta" = "#AB47BC", "purple" = "#7B1FA2",
               "grey" = "#9E9E9E")

  pb <- ggplot(kda_top_per_mod,
               aes(x = -log10(padj), y = gene, fill = module)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("k=%d, FE=%.1f", module_neighbors, fold_enrichment)),
              hjust = -0.05, size = 2.3) +
    geom_point(data = kda_top_per_mod[is_deg == TRUE],
               aes(x = -log10(padj) * 1.02), shape = 18, size = 2,
               color = masld_colors$up, show.legend = FALSE) +
    geom_point(data = kda_top_per_mod[is_coloc == TRUE],
               aes(x = -log10(padj) * 1.04), shape = 17, size = 2,
               color = masld_colors$gwas, show.legend = FALSE) +
    scale_fill_manual(values = mod_pal, name = "WGCNA\nModule") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(x = expression(-log[10](p[adj])),
         y = NULL,
         title = "Top Key Drivers per WGCNA Module") +
    theme_masld() +
    theme(axis.text.y = element_text(size = 7))
} else {
  pb <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No significant key drivers\n(padj < 0.05)", size = 5) +
    theme_void() +
    ggtitle("Key Driver Analysis")
}

# --- Panel C: Network of top connectors linking DEG + COLOC seeds ---
# Extract subnetwork of DIAMOnD connector genes + their DEG/COLOC neighbors
top_connectors <- diamond_results[category == "Connector"][order(p_value)][1:min(30, .N)]$gene
# Include top COLOC genes in PPI for the visualization
seed_nodes <- unique(c(
  head(top500_in_ppi, 30),          # top DEG seeds
  head(coloc_in_ppi, 20),           # COLOC seeds
  top_connectors                     # top connector genes
))
seed_nodes <- intersect(seed_nodes, all_genes_ppi)

if (length(seed_nodes) > 5) {
  sub_g <- induced_subgraph(g, seed_nodes)
  # Remove isolates
  sub_g <- delete_vertices(sub_g, which(degree(sub_g) == 0))

  if (vcount(sub_g) > 2) {
    # Classify nodes
    node_type <- ifelse(
      V(sub_g)$name %in% deg_symbols & V(sub_g)$name %in% coloc_genes, "DEG+COLOC",
      ifelse(V(sub_g)$name %in% deg_symbols, "DEG",
      ifelse(V(sub_g)$name %in% coloc_genes, "COLOC", "Connector"))
    )

    # Layout
    set.seed(42)
    lo <- layout_with_fr(sub_g)
    net_df <- data.table(
      x = lo[, 1], y = lo[, 2],
      name = V(sub_g)$name,
      type = node_type,
      degree = degree(sub_g)
    )

    edge_list <- as_edgelist(sub_g)
    edge_df <- data.table(
      x = lo[match(edge_list[, 1], V(sub_g)$name), 1],
      y = lo[match(edge_list[, 1], V(sub_g)$name), 2],
      xend = lo[match(edge_list[, 2], V(sub_g)$name), 1],
      yend = lo[match(edge_list[, 2], V(sub_g)$name), 2]
    )

    type_colors <- c(
      "DEG" = masld_colors$up,
      "COLOC" = masld_colors$gwas,
      "DEG+COLOC" = masld_colors$fibrosis,
      "Connector" = "#66BB6A"
    )

    pc <- ggplot() +
      geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
                   color = "grey80", linewidth = 0.3, alpha = 0.6) +
      geom_point(data = net_df, aes(x = x, y = y, color = type, size = degree)) +
      geom_text(data = net_df[degree >= quantile(degree, 0.6)],
                aes(x = x, y = y, label = name),
                size = 2.2, vjust = -1, check_overlap = TRUE) +
      scale_color_manual(values = type_colors, name = "Gene Type") +
      scale_size_continuous(range = c(1.5, 5), name = "Degree") +
      labs(title = "Connector Network: DEG + COLOC Seeds") +
      theme_masld() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank())
  } else {
    pc <- ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "Insufficient connected nodes", size = 5) +
      theme_void()
  }
} else {
  pc <- ggplot() + annotate("text", x = 0.5, y = 0.5,
                             label = "Insufficient seed nodes in PPI", size = 5) +
    theme_void()
}

# --- Compose ---
fig <- pa_final / (pb | pc) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )

ggsave(out_fig, fig, width = 14, height = 12, units = "in")
cat(sprintf("  Figure saved: %s\n", out_fig))

cat("\nDone:", format(Sys.time()), "\n")
