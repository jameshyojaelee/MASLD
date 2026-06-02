#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_256_spatial_cosmos_xsp
#SBATCH --output=logs/net_256_spatial_cosmos_xsp_%j.out
#SBATCH --error=logs/net_256_spatial_cosmos_xsp_%j.err
# ===========================================================================
# Script 256: Build spatial, COSMOS, and cross-species edges for the
#             Bayesian multiplex gene network (Layers 7, 8, 10).
# ===========================================================================
#
# Input:
#   - RNA-seq/results/network/network_nodes.csv  (node set V from Script 250)
#   - Analysis/Spatial/results/coexpression/spatial_autocorr_{Healthy,Steatotic}.csv
#   - Analysis/Spatial/results/coexpression/modules_{Healthy,Steatotic}.csv
#   - Analysis/Spatial/results/svg/differential_svgs.csv
#   - RNA-seq/results/multi_evidence/cosmos_mechanistic/gwas_tf_de_paths.csv
#   - RNA-seq/results/multi_evidence/cosmos_mechanistic/tf_gwas_de_bridge.csv
#   - Analysis/Cross_Species_Concordance/results/gene_concordance_per_gene.csv
#   - RNA-seq/Mouse/Unified_Integration/results/per_diet/{diet}_de_results.csv
#   - streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz
#
# Output:
#   - RNA-seq/results/network/edges_spatial.parquet   (Layer 7)
#   - RNA-seq/results/network/edges_cosmos.parquet    (Layer 8)
#   - RNA-seq/results/network/edges_xspecies.parquet  (Layer 10)
#
# Environment: rnaseq (R 4.4+, data.table, arrow)
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  # library(arrow)  # using fwrite instead
})

t0 <- Sys.time()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR <- file.path(BASE, "RNA-seq/results/network")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

NODE_PATH <- file.path(OUTDIR, "network_nodes.csv")
stopifnot(file.exists(NODE_PATH))

cat("=== 256: Spatial / COSMOS / Cross-Species Edge Construction ===\n")
cat("Start:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---------------------------------------------------------------------------
# 0. Load node set V
# ---------------------------------------------------------------------------
nodes <- fread(NODE_PATH)
V <- unique(nodes$human_symbol)
V <- V[!is.na(V) & V != ""]
cat("Node set V:", length(V), "genes\n\n")

# ===========================================================================
#  LAYER 7 — Spatial co-localization edges
# ===========================================================================
cat("###########################################################\n")
cat("# LAYER 7: Spatial Co-localization (GSE192741 Visium)     #\n")
cat("###########################################################\n\n")

t7 <- Sys.time()

SPATIAL_DIR <- file.path(BASE, "Analysis/Spatial/results")

# --- 7a. Load Moran's I spatial autocorrelation for both conditions ---
# Files have columns: Gene, C (Moran's C), Z, Pval, FDR
autocorr_h_path <- file.path(SPATIAL_DIR, "coexpression/spatial_autocorr_Healthy.csv")
autocorr_s_path <- file.path(SPATIAL_DIR, "coexpression/spatial_autocorr_Steatotic.csv")
mod_h_path      <- file.path(SPATIAL_DIR, "coexpression/modules_Healthy.csv")
mod_s_path      <- file.path(SPATIAL_DIR, "coexpression/modules_Steatotic.csv")
svg_path        <- file.path(SPATIAL_DIR, "svg/differential_svgs.csv")

stopifnot(file.exists(autocorr_h_path), file.exists(autocorr_s_path))
stopifnot(file.exists(mod_h_path), file.exists(mod_s_path))
stopifnot(file.exists(svg_path))

# Load spatial autocorrelation scores (Moran's C ~ spatial co-variation)
autocorr_h <- fread(autocorr_h_path)
autocorr_s <- fread(autocorr_s_path)
cat("Spatial autocorrelation: Healthy", nrow(autocorr_h), "genes, Steatotic", nrow(autocorr_s), "genes\n")

# Load module assignments (co-expression modules from spatial data)
# Columns: V1 (gene symbol), Module
mod_h <- fread(mod_h_path)
mod_s <- fread(mod_s_path)
setnames(mod_h, c("gene", "module_healthy"))
setnames(mod_s, c("gene", "module_steatotic"))
cat("Spatial modules: Healthy", nrow(mod_h), "genes (", uniqueN(mod_h$module_healthy), "modules),",
    "Steatotic", nrow(mod_s), "genes (", uniqueN(mod_s$module_steatotic), "modules)\n")

# Load differential SVGs (genes with Moran's I in both conditions)
# Columns: V1 (gene), morans_I_healthy, morans_I_masld, svg_healthy, svg_masld, delta_I, category
svgs <- fread(svg_path)
setnames(svgs, "V1", "gene")
cat("Differential SVGs:", nrow(svgs), "genes\n")

# --- 7b. Build spatial profile for each gene ---
# Strategy: For each gene, build a profile vector:
#   [Moran's_C_healthy, Moran's_C_steatotic, morans_I_healthy, morans_I_masld,
#    module_healthy, module_steatotic]
# Then compute pairwise similarity between genes in V.

# Merge autocorrelation data
setnames(autocorr_h, c("Gene", "C_healthy", "Z_healthy", "Pval_healthy", "FDR_healthy"))
setnames(autocorr_s, c("Gene", "C_steatotic", "Z_steatotic", "Pval_steatotic", "FDR_steatotic"))

spatial_prof <- merge(autocorr_h[, .(Gene, C_healthy)],
                      autocorr_s[, .(Gene, C_steatotic)],
                      by = "Gene", all = FALSE)

# Merge Moran's I from SVG differential file
spatial_prof <- merge(spatial_prof,
                      svgs[, .(gene, morans_I_healthy, morans_I_masld)],
                      by.x = "Gene", by.y = "gene", all.x = TRUE)

# Merge module assignments
spatial_prof <- merge(spatial_prof, mod_h, by.x = "Gene", by.y = "gene", all.x = TRUE)
spatial_prof <- merge(spatial_prof, mod_s, by.x = "Gene", by.y = "gene", all.x = TRUE)

# Filter to genes in V
spatial_prof <- spatial_prof[Gene %in% V]
cat("Genes in V with spatial data:", nrow(spatial_prof), "\n")

if (nrow(spatial_prof) >= 2) {

  # --- 7c. Compute pairwise Pearson on continuous spatial features ---
  # Use Moran's C (spatial autocorrelation coefficient) from both conditions
  # as the spatial expression profile.
  feat_cols <- c("C_healthy", "C_steatotic")
  if ("morans_I_healthy" %in% names(spatial_prof)) {
    # Add Moran's I if available (from SVG analysis)
    feat_cols <- c(feat_cols, "morans_I_healthy", "morans_I_masld")
  }

  # Keep only genes with complete cases on continuous features
  complete_mask <- complete.cases(spatial_prof[, ..feat_cols])
  sp_complete <- spatial_prof[complete_mask]
  cat("Genes with complete spatial features:", nrow(sp_complete), "\n")

  if (nrow(sp_complete) >= 2) {

    # Build numeric matrix (genes x features)
    sp_mat <- as.matrix(sp_complete[, ..feat_cols])
    rownames(sp_mat) <- sp_complete$Gene

    n_sp <- nrow(sp_mat)
    cat("Computing", n_sp, "x", n_sp, "Pearson correlation on spatial profiles ...\n")

    cor_mat <- cor(t(sp_mat), method = "pearson", use = "pairwise.complete.obs")

    # Zero out lower triangle + diagonal
    cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA

    # Determine top 1% threshold by |r|
    upper_vals <- cor_mat[upper.tri(cor_mat, diag = FALSE)]
    upper_vals <- upper_vals[!is.na(upper_vals)]
    n_pairs <- length(upper_vals)
    thresh_spatial <- quantile(abs(upper_vals), probs = 0.99, na.rm = TRUE)
    cat("Total unique pairs:", format(n_pairs, big.mark = ","), "\n")
    cat("|r| threshold at top 1%:", round(thresh_spatial, 4), "\n")

    rm(upper_vals); gc(verbose = FALSE)

    # Extract edges above threshold
    idx <- which(abs(cor_mat) >= thresh_spatial, arr.ind = TRUE)
    gene_names_sp <- rownames(sp_mat)

    edges_sp <- data.table(
      gene_a    = gene_names_sp[idx[, 1]],
      gene_b    = gene_names_sp[idx[, 2]],
      raw_score = cor_mat[idx]
    )
    rm(cor_mat, idx); gc(verbose = FALSE)

    # Remove self-loops, canonical ordering
    edges_sp <- edges_sp[gene_a != gene_b]
    swap <- edges_sp$gene_a > edges_sp$gene_b
    if (any(swap)) {
      tmp_a <- edges_sp$gene_a[swap]
      edges_sp[swap, gene_a := gene_b]
      edges_sp[swap, gene_b := tmp_a]
      rm(tmp_a)
    }
    edges_sp <- unique(edges_sp, by = c("gene_a", "gene_b"))

    # --- 7d. Also add co-module edges (same spatial module = similar spatial niche) ---
    # Genes in the same spatial co-expression module share spatial patterns.
    # Create edges for gene pairs sharing a module in EITHER condition.
    cat("\n--- Adding spatial co-module edges ---\n")

    make_comod_edges <- function(mod_dt, mod_col, condition_label) {
      mod_v <- mod_dt[gene %in% V & !is.na(get(mod_col))]
      if (nrow(mod_v) < 2) return(data.table())

      # For each module, create edges between its member genes
      mod_list <- split(mod_v$gene, mod_v[[mod_col]])
      edge_list <- lapply(mod_list, function(genes) {
        if (length(genes) < 2) return(NULL)
        # Use combn for pairwise combinations
        pairs <- combn(sort(genes), 2)
        data.table(gene_a = pairs[1, ], gene_b = pairs[2, ],
                   module = condition_label)
      })
      rbindlist(edge_list[!sapply(edge_list, is.null)])
    }

    comod_h <- make_comod_edges(mod_h, "module_healthy", "healthy")
    comod_s <- make_comod_edges(mod_s, "module_steatotic", "steatotic")
    cat("Co-module edges: Healthy", nrow(comod_h), ", Steatotic", nrow(comod_s), "\n")

    # Merge co-module edges: if a pair shares a module in both conditions, note it
    comod_all <- rbind(comod_h, comod_s)
    if (nrow(comod_all) > 0) {
      comod_summary <- comod_all[, .(
        n_shared_modules = .N,
        conditions = paste(unique(module), collapse = ";")
      ), by = .(gene_a, gene_b)]

      # Assign score: 1.0 if shared in both conditions, 0.5 if one condition
      comod_summary[, comod_score := fifelse(n_shared_modules >= 2, 1.0, 0.5)]
      cat("Unique co-module pairs:", nrow(comod_summary), "\n")
    } else {
      comod_summary <- data.table(gene_a = character(), gene_b = character(),
                                  comod_score = numeric(), conditions = character())
    }

    # --- 7e. Merge correlation-based and co-module edges ---
    # Correlation edges get their r as raw_score; co-module edges get comod_score
    # For genes with both types, take max score

    if (nrow(edges_sp) > 0) {
      # Compute p-values for correlation edges (t-distribution approximation)
      # Here n = number of features (small), so p-values are approximate
      n_feat <- length(feat_cols)
      df_sp <- n_feat - 2L
      if (df_sp >= 1) {
        edges_sp[, t_stat := raw_score * sqrt(df_sp / (1 - raw_score^2 + 1e-10))]
        edges_sp[, pvalue := 2 * pt(-abs(t_stat), df = df_sp)]
      } else {
        edges_sp[, pvalue := NA_real_]
      }

      edges_sp[, metadata := paste0("p=", signif(pvalue, 4), ";dataset=GSE192741")]
      edges_sp[, source := "correlation"]
    }

    if (nrow(comod_summary) > 0) {
      comod_edges <- comod_summary[, .(
        gene_a, gene_b,
        raw_score = comod_score,
        metadata = paste0("p=NA;dataset=GSE192741;modules=", conditions),
        source = "comodule"
      )]
    } else {
      comod_edges <- data.table(gene_a = character(), gene_b = character(),
                                raw_score = numeric(), metadata = character(),
                                source = character())
    }

    # Combine, keeping max score per pair
    all_spatial <- rbind(
      edges_sp[, .(gene_a, gene_b, raw_score, metadata, source)],
      comod_edges,
      fill = TRUE
    )
    all_spatial[, abs_score := abs(raw_score)]
    setorder(all_spatial, gene_a, gene_b, -abs_score)
    all_spatial <- all_spatial[!duplicated(all_spatial[, .(gene_a, gene_b)])]
    all_spatial[, abs_score := NULL]
    all_spatial[, source := NULL]

  } else {
    cat("WARN: Too few genes with complete spatial features for correlation\n")
    all_spatial <- data.table(gene_a = character(), gene_b = character(),
                              raw_score = numeric(), metadata = character())
  }
} else {
  cat("WARN: Too few genes in V with spatial data\n")
  all_spatial <- data.table(gene_a = character(), gene_b = character(),
                            raw_score = numeric(), metadata = character())
}

# Save Layer 7
out_sp <- all_spatial[, .(gene_a, gene_b, raw_score, metadata)]
setorder(out_sp, -abs(raw_score))
fwrite(out_sp, file.path(OUTDIR, "edges_spatial.csv"))

cat("\n--- Layer 7 Summary ---\n")
cat("Total spatial edges:", format(nrow(out_sp), big.mark = ","), "\n")
if (nrow(out_sp) > 0) {
  cat("raw_score range: [", round(min(out_sp$raw_score), 4), ",",
      round(max(out_sp$raw_score), 4), "]\n")
  cat("Mean |raw_score|:", round(mean(abs(out_sp$raw_score)), 4), "\n")
  cat("Median |raw_score|:", round(median(abs(out_sp$raw_score)), 4), "\n")

  # Top genes by degree
  deg_a <- out_sp[, .N, by = gene_a]; setnames(deg_a, c("gene", "n"))
  deg_b <- out_sp[, .N, by = gene_b]; setnames(deg_b, c("gene", "n"))
  deg_all <- rbind(deg_a, deg_b)[, .(degree = sum(n)), by = gene]
  setorder(deg_all, -degree)
  cat("Top 10 genes by spatial degree:\n")
  print(head(deg_all, 10))
}
cat("Layer 7 time:", round(difftime(Sys.time(), t7, units = "secs"), 1), "s\n\n")

# ===========================================================================
#  LAYER 8 — COSMOS mechanistic path edges
# ===========================================================================
cat("###########################################################\n")
cat("# LAYER 8: COSMOS Mechanistic Paths                       #\n")
cat("###########################################################\n\n")

t8 <- Sys.time()

COSMOS_DIR <- file.path(BASE, "RNA-seq/results/multi_evidence/cosmos_mechanistic")
paths_path  <- file.path(COSMOS_DIR, "gwas_tf_de_paths.csv")
bridge_path <- file.path(COSMOS_DIR, "tf_gwas_de_bridge.csv")

stopifnot(file.exists(paths_path), file.exists(bridge_path))

# --- 8a. Load COSMOS paths ---
# gwas_tf_de_paths.csv columns:
#   tf, tf_activity, gwas_gene (semicolon-delimited), n_de_targets,
#   top_de_targets (semicolon-delimited), path_description
cosmos_paths <- fread(paths_path)
cat("COSMOS paths:", nrow(cosmos_paths), "TFs\n")

# tf_gwas_de_bridge.csv columns:
#   tf, tf_activity, n_total_targets, n_gwas_targets, expected_gwas,
#   gwas_enrichment, gwas_hypergeom_pval, n_deg_targets, gwas_targets (;-sep),
#   bridges_gwas_to_de, gwas_hypergeom_padj
cosmos_bridge <- fread(bridge_path)
cat("COSMOS bridge TFs:", nrow(cosmos_bridge), "\n")

# --- 8b. Build edges from paths ---
# Three edge types from COSMOS:
# 1. gwas_to_tf: GWAS variant gene -> TF (path_length = 2, score = 0.5)
# 2. tf_to_de: TF -> downstream DE target (path_length = 2, score = 0.5)
# 3. direct: GWAS gene -> DE target via TF (path_length = 3, score = 0.33)

edge_list_cosmos <- list()

for (i in seq_len(nrow(cosmos_paths))) {
  tf <- cosmos_paths$tf[i]
  gwas_genes <- unlist(strsplit(cosmos_paths$gwas_gene[i], ";"))
  gwas_genes <- trimws(gwas_genes)
  gwas_genes <- gwas_genes[gwas_genes != ""]

  de_targets <- unlist(strsplit(cosmos_paths$top_de_targets[i], ";"))
  de_targets <- trimws(de_targets)
  de_targets <- de_targets[de_targets != ""]

  # Edge type 1: GWAS gene -> TF (path_length = 2)
  gwas_in_v <- gwas_genes[gwas_genes %in% V]
  if (tf %in% V && length(gwas_in_v) > 0) {
    for (g in gwas_in_v) {
      pair <- sort(c(g, tf))
      edge_list_cosmos[[length(edge_list_cosmos) + 1]] <- data.table(
        gene_a = pair[1], gene_b = pair[2],
        raw_score = 0.5,
        path_type = "gwas_to_tf",
        tf_node = tf
      )
    }
  }

  # Edge type 2: TF -> DE target (path_length = 2)
  de_in_v <- de_targets[de_targets %in% V]
  if (tf %in% V && length(de_in_v) > 0) {
    for (d in de_in_v) {
      pair <- sort(c(tf, d))
      edge_list_cosmos[[length(edge_list_cosmos) + 1]] <- data.table(
        gene_a = pair[1], gene_b = pair[2],
        raw_score = 0.5,
        path_type = "tf_to_de",
        tf_node = tf
      )
    }
  }

  # Edge type 3: Direct GWAS gene -> DE target (path_length = 3, through TF)
  if (length(gwas_in_v) > 0 && length(de_in_v) > 0) {
    for (g in gwas_in_v) {
      for (d in de_in_v) {
        if (g == d) next
        pair <- sort(c(g, d))
        edge_list_cosmos[[length(edge_list_cosmos) + 1]] <- data.table(
          gene_a = pair[1], gene_b = pair[2],
          raw_score = 1 / 3,
          path_type = "direct",
          tf_node = tf
        )
      }
    }
  }
}

if (length(edge_list_cosmos) > 0) {
  edges_cosmos <- rbindlist(edge_list_cosmos)
  cat("Raw COSMOS edges (before dedup):", format(nrow(edges_cosmos), big.mark = ","), "\n")

  # Remove self-loops
  edges_cosmos <- edges_cosmos[gene_a != gene_b]

  # For duplicate pairs, keep the highest raw_score (shortest path) and
  # concatenate path_types
  edges_cosmos_dedup <- edges_cosmos[, .(
    raw_score = max(raw_score),
    path_type = paste(unique(path_type), collapse = ";"),
    tf_nodes  = paste(unique(tf_node), collapse = ";")
  ), by = .(gene_a, gene_b)]

  # Build metadata
  edges_cosmos_dedup[, metadata := paste0("path_type=", path_type, ";tf=", tf_nodes)]
  edges_cosmos_dedup[, c("path_type", "tf_nodes") := NULL]
} else {
  edges_cosmos_dedup <- data.table(gene_a = character(), gene_b = character(),
                                   raw_score = numeric(), metadata = character())
}

# Save Layer 8
out_cosmos <- edges_cosmos_dedup[, .(gene_a, gene_b, raw_score, metadata)]
setorder(out_cosmos, -raw_score)
fwrite(out_cosmos, file.path(OUTDIR, "edges_cosmos.csv"))

cat("\n--- Layer 8 Summary ---\n")
cat("Total COSMOS edges:", format(nrow(out_cosmos), big.mark = ","), "\n")
if (nrow(out_cosmos) > 0) {
  cat("raw_score range: [", round(min(out_cosmos$raw_score), 4), ",",
      round(max(out_cosmos$raw_score), 4), "]\n")
  cat("Mean raw_score:", round(mean(out_cosmos$raw_score), 4), "\n")

  # Path type breakdown
  all_types <- unlist(strsplit(out_cosmos$metadata, ";"))
  type_vals <- grep("^path_type=", all_types, value = TRUE)
  type_vals <- sub("^path_type=", "", type_vals)
  type_counts <- table(unlist(strsplit(type_vals, ";")))
  cat("Path type breakdown:\n")
  print(type_counts)

  # Top genes by degree
  deg_a <- out_cosmos[, .N, by = gene_a]; setnames(deg_a, c("gene", "n"))
  deg_b <- out_cosmos[, .N, by = gene_b]; setnames(deg_b, c("gene", "n"))
  deg_all <- rbind(deg_a, deg_b)[, .(degree = sum(n)), by = gene]
  setorder(deg_all, -degree)
  cat("\nTop 10 genes by COSMOS degree:\n")
  print(head(deg_all, 10))
}
cat("Layer 8 time:", round(difftime(Sys.time(), t8, units = "secs"), 1), "s\n\n")

# ===========================================================================
#  LAYER 10 — Cross-species concordance edges
# ===========================================================================
cat("###########################################################\n")
cat("# LAYER 10: Cross-Species Concordance                     #\n")
cat("###########################################################\n\n")

t10 <- Sys.time()

XSPEC_DIR <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
conc_path <- file.path(XSPEC_DIR, "gene_concordance_per_gene.csv")
MOUSE_DIET_DIR <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
ORTHO_PATH <- file.path(BASE, "streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz")

stopifnot(file.exists(conc_path))

# --- 10a. Load concordance data ---
# Columns: mouse_gene_id, human_symbol, h_significant, n_mouse_sig, n_diets_sig,
#   n_concordant, n_discordant, mean_h_lfc, diets_concordant, diets_discordant,
#   diets_mouse_sig, primary_category, best_n_concordant, best_category, ...
conc <- fread(conc_path)
cat("Concordance data:", nrow(conc), "gene pairs\n")
cat("Categories:", paste(sort(unique(conc$primary_category)), collapse = ", "), "\n")

# Filter to genes in V
conc_v <- conc[human_symbol %in% V]
cat("Genes in V with concordance data:", nrow(conc_v), "\n")

# --- 10b. Build per-diet logFC matrix for concordance similarity ---
# Load per-diet DE results to get diet-specific logFC for each gene.
# This lets us compute correlation of logFC profiles across diets.
DIETS <- c("MCD", "HFD", "CDAHFD", "FPC", "LIDPAD")
diet_lfc_list <- list()

for (diet in DIETS) {
  f <- file.path(MOUSE_DIET_DIR, paste0(diet, "_de_results.csv"))
  if (!file.exists(f)) {
    cat("  WARN: Missing diet file:", f, "\n")
    next
  }
  dt <- fread(f, select = c("gene", "logFC"))
  setnames(dt, c("mouse_gene_id", paste0("lfc_", diet)))
  diet_lfc_list[[diet]] <- dt
}

if (length(diet_lfc_list) > 0) {
  # Merge all diets into one matrix
  diet_lfc <- Reduce(function(x, y) merge(x, y, by = "mouse_gene_id", all = TRUE),
                     diet_lfc_list)
  cat("Mouse diet logFC matrix:", nrow(diet_lfc), "genes x",
      ncol(diet_lfc) - 1, "diets\n")

  # Merge with concordance data to get human symbols
  conc_v_lfc <- merge(conc_v[, .(mouse_gene_id, human_symbol, primary_category,
                                  n_concordant, n_discordant, diets_concordant,
                                  n_diets_sig)],
                      diet_lfc,
                      by = "mouse_gene_id", all.x = TRUE)
  cat("Concordance genes with diet logFC:", sum(complete.cases(
    conc_v_lfc[, .SD, .SDcols = grep("^lfc_", names(conc_v_lfc), value = TRUE)]
  )), "\n")
} else {
  cat("WARN: No diet logFC files found; falling back to category-only similarity\n")
  conc_v_lfc <- conc_v[, .(mouse_gene_id, human_symbol, primary_category,
                            n_concordant, n_discordant, diets_concordant,
                            n_diets_sig)]
}

# --- 10c. Create edges for genes sharing concordance category ---
# Concordant categories that indicate shared biology:
# Conserved, High_Concordance, Moderate_Concordance, Species_Discordant
# Edges connect genes in V that share the same non-trivial category.

active_categories <- c("Conserved", "High_Concordance", "Moderate_Concordance")

edge_list_xsp <- list()

for (cat_name in active_categories) {
  genes_in_cat <- conc_v_lfc[primary_category == cat_name, human_symbol]
  genes_in_cat <- unique(genes_in_cat[genes_in_cat %in% V])

  if (length(genes_in_cat) < 2) next

  # Sort for canonical ordering
  genes_in_cat <- sort(genes_in_cat)
  n_genes <- length(genes_in_cat)
  cat("  Category", cat_name, ":", n_genes, "genes in V ->",
      choose(n_genes, 2), "potential pairs\n")

  # Cap pairwise edges for very large categories to avoid memory blow-up
  MAX_GENES_PAIRWISE <- 2000
  if (n_genes > MAX_GENES_PAIRWISE) {
    cat("    Capping to top", MAX_GENES_PAIRWISE,
        "genes by n_concordant for pairwise edges\n")
    top_genes <- conc_v_lfc[primary_category == cat_name & human_symbol %in% V]
    setorder(top_genes, -n_concordant)
    genes_in_cat <- sort(unique(top_genes$human_symbol[seq_len(MAX_GENES_PAIRWISE)]))
    n_genes <- length(genes_in_cat)
    cat("    Adjusted:", n_genes, "genes ->", choose(n_genes, 2), "potential pairs\n")
  }

  # Compute pairwise similarity using diet logFC correlation
  lfc_cols <- grep("^lfc_", names(conc_v_lfc), value = TRUE)

  if (length(lfc_cols) > 0) {
    # Build logFC matrix for genes in this category
    lfc_sub <- conc_v_lfc[human_symbol %in% genes_in_cat & !duplicated(human_symbol)]
    lfc_mat <- as.matrix(lfc_sub[, ..lfc_cols])
    rownames(lfc_mat) <- lfc_sub$human_symbol

    # Remove genes with all NA logFC
    complete_mask <- rowSums(!is.na(lfc_mat)) >= 2
    lfc_mat <- lfc_mat[complete_mask, , drop = FALSE]

    if (nrow(lfc_mat) >= 2) {
      # Compute Pearson correlation of logFC profiles across diets
      cor_xsp <- cor(t(lfc_mat), method = "pearson", use = "pairwise.complete.obs")
      cor_xsp[lower.tri(cor_xsp, diag = TRUE)] <- NA

      # Extract all pairs with positive correlation (same direction in mouse)
      idx <- which(!is.na(cor_xsp) & cor_xsp > 0, arr.ind = TRUE)
      gene_names_xsp <- rownames(lfc_mat)

      if (nrow(idx) > 0) {
        pairs_dt <- data.table(
          gene_a    = gene_names_xsp[idx[, 1]],
          gene_b    = gene_names_xsp[idx[, 2]],
          raw_score = cor_xsp[idx],
          concordance_category = cat_name
        )
        edge_list_xsp[[length(edge_list_xsp) + 1]] <- pairs_dt
      }
      rm(cor_xsp, idx); gc(verbose = FALSE)
    }
  } else {
    # Fallback: If no logFC matrix, create edges with score based on shared
    # concordant diet count (Jaccard-like)
    for (ii in seq_len(n_genes - 1)) {
      gene_i <- genes_in_cat[ii]
      diets_i <- unlist(strsplit(
        conc_v_lfc[human_symbol == gene_i, diets_concordant][1], ";"
      ))
      diets_i <- diets_i[diets_i != ""]

      for (jj in (ii + 1):n_genes) {
        gene_j <- genes_in_cat[jj]
        diets_j <- unlist(strsplit(
          conc_v_lfc[human_symbol == gene_j, diets_concordant][1], ";"
        ))
        diets_j <- diets_j[diets_j != ""]

        if (length(diets_i) == 0 || length(diets_j) == 0) next

        shared <- length(intersect(diets_i, diets_j))
        union_n <- length(union(diets_i, diets_j))
        jaccard <- shared / union_n

        if (jaccard > 0) {
          edge_list_xsp[[length(edge_list_xsp) + 1]] <- data.table(
            gene_a = gene_i, gene_b = gene_j,
            raw_score = jaccard,
            concordance_category = cat_name
          )
        }
      }
    }
  }
}

if (length(edge_list_xsp) > 0) {
  edges_xsp <- rbindlist(edge_list_xsp)
  cat("\nRaw cross-species edges:", format(nrow(edges_xsp), big.mark = ","), "\n")

  # Deduplicate: for pairs appearing in multiple categories, keep max score
  edges_xsp <- edges_xsp[gene_a != gene_b]
  swap <- edges_xsp$gene_a > edges_xsp$gene_b
  if (any(swap)) {
    tmp_a <- edges_xsp$gene_a[swap]
    edges_xsp[swap, gene_a := gene_b]
    edges_xsp[swap, gene_b := tmp_a]
    rm(tmp_a)
  }

  edges_xsp_dedup <- edges_xsp[, .(
    raw_score = max(raw_score),
    concordance_categories = paste(unique(concordance_category), collapse = ";")
  ), by = .(gene_a, gene_b)]

  # Build shared_diets metadata by looking up each gene's concordant diets
  diet_lookup <- conc_v_lfc[, .(
    human_symbol, diets_concordant, n_diets_sig
  )][!duplicated(human_symbol)]

  edges_xsp_dedup <- merge(edges_xsp_dedup,
                           diet_lookup[, .(human_symbol, diets_a = diets_concordant)],
                           by.x = "gene_a", by.y = "human_symbol", all.x = TRUE)
  edges_xsp_dedup <- merge(edges_xsp_dedup,
                           diet_lookup[, .(human_symbol, diets_b = diets_concordant)],
                           by.x = "gene_b", by.y = "human_symbol", all.x = TRUE)

  edges_xsp_dedup[, shared_diets := mapply(function(a, b) {
    da <- unlist(strsplit(as.character(a), ";"))
    db <- unlist(strsplit(as.character(b), ";"))
    da <- da[da != "" & !is.na(da)]
    db <- db[db != "" & !is.na(db)]
    paste(intersect(da, db), collapse = ";")
  }, diets_a, diets_b)]

  edges_xsp_dedup[, metadata := paste0(
    "shared_diets=", fifelse(shared_diets == "", "none", shared_diets),
    ";concordance_categories=", concordance_categories
  )]
  edges_xsp_dedup[, c("diets_a", "diets_b", "shared_diets",
                       "concordance_categories") := NULL]
} else {
  edges_xsp_dedup <- data.table(gene_a = character(), gene_b = character(),
                                raw_score = numeric(), metadata = character())
}

# Save Layer 10
out_xsp <- edges_xsp_dedup[, .(gene_a, gene_b, raw_score, metadata)]
setorder(out_xsp, -raw_score)
fwrite(out_xsp, file.path(OUTDIR, "edges_xspecies.csv"))

cat("\n--- Layer 10 Summary ---\n")
cat("Total cross-species edges:", format(nrow(out_xsp), big.mark = ","), "\n")
if (nrow(out_xsp) > 0) {
  cat("raw_score range: [", round(min(out_xsp$raw_score), 4), ",",
      round(max(out_xsp$raw_score), 4), "]\n")
  cat("Mean raw_score:", round(mean(out_xsp$raw_score), 4), "\n")
  cat("Median raw_score:", round(median(out_xsp$raw_score), 4), "\n")

  # Top genes by degree
  deg_a <- out_xsp[, .N, by = gene_a]; setnames(deg_a, c("gene", "n"))
  deg_b <- out_xsp[, .N, by = gene_b]; setnames(deg_b, c("gene", "n"))
  deg_all <- rbind(deg_a, deg_b)[, .(degree = sum(n)), by = gene]
  setorder(deg_all, -degree)
  cat("\nTop 10 genes by cross-species degree:\n")
  print(head(deg_all, 10))

  # Category breakdown
  cat_vals <- gsub(".*concordance_categories=", "", out_xsp$metadata)
  cat("\nEdges by concordance category:\n")
  print(table(cat_vals))
}
cat("Layer 10 time:", round(difftime(Sys.time(), t10, units = "secs"), 1), "s\n\n")

# ===========================================================================
# Final summary
# ===========================================================================
cat("###########################################################\n")
cat("# FINAL SUMMARY                                           #\n")
cat("###########################################################\n\n")

cat(sprintf("%-25s %10s  %s\n", "Layer", "Edges", "Output"))
cat(strrep("-", 65), "\n")
cat(sprintf("%-25s %10s  %s\n", "L7: Spatial",
            format(nrow(out_sp), big.mark = ","),
            "edges_spatial.parquet"))
cat(sprintf("%-25s %10s  %s\n", "L8: COSMOS",
            format(nrow(out_cosmos), big.mark = ","),
            "edges_cosmos.parquet"))
cat(sprintf("%-25s %10s  %s\n", "L10: Cross-species",
            format(nrow(out_xsp), big.mark = ","),
            "edges_xspecies.parquet"))

t1 <- Sys.time()
cat("\nTotal elapsed:", round(difftime(t1, t0, units = "secs"), 1), "s\n")
cat("Done.\n")
