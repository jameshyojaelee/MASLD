#!/usr/bin/env Rscript
# ===========================================================================
# Script 59: Multiplex Network Random Walk with Restart (RWR)
# ===========================================================================
# Purpose: Score ALL genes (including lncRNAs) by network proximity to
#          high-confidence seed genes (DEG + COLOC overlap) using RWR
#          across a multiplex network (PPI + coexpression + regulatory).
#
# Layers:
#   L1 — STRING PPI (combined_score > 700)
#   L2 — Coexpression (top 0.1% pairwise Pearson from merged logCPM)
#   L3 — Regulatory (WGCNA module co-membership edges)
#
# Seeds: 196 DEG+COLOC overlap genes (dream padj<0.05, |logFC|>0.3,
#         COLOC PP.H4>0.5)
#
# RWR: p(t+1) = (1-r)*M*p(t) + r*p(0), r=0.7, convergence tol=1e-8
#
# Outputs:
#   - RNA-seq/results/multi_evidence/rwr_scores.csv
#   - figures/supplementary/figS08_subtyping_convergence/figS_multiplex_rwr.pdf
#
# Environment: rnaseq (R 4.4+, igraph, Matrix, data.table, ggplot2)
# SLURM: cpu partition, 8 CPUs, 64G, 48h
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(igraph)
  library(ggplot2)
  library(scales)
  library(parallel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

cat("=== Script 59: Multiplex RWR ===\n")
cat("Start:", format(Sys.time()), "\n")

# ---------------------------------------------------------------------------
# 0. Output paths
# ---------------------------------------------------------------------------
out_csv <- file.path(BASE, "RNA-seq/results/multi_evidence/rwr_scores.csv")
fig_dir <- file.path(BASE, "figures/supplementary/figS08_subtyping_convergence")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(fig_dir, "figS_multiplex_rwr.pdf")

# ---------------------------------------------------------------------------
# 1. Load multi-evidence atlas for seeds + gene metadata
# ---------------------------------------------------------------------------
cat("\n--- Loading atlas ---\n")
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
cat("Atlas:", nrow(atlas), "genes x", ncol(atlas), "columns\n")

# Identify seed genes: DEG (padj<0.05, |logFC|>0.3) AND COLOC PP.H4>0.5
# Use coloc_susie_best_pp4 (SuSiE-COLOC, active column with 577 genes > 0.5)
atlas[, is_deg := !is.na(dream_padj) & dream_padj < 0.05 & abs(dream_logFC) > 0.3]
atlas[, has_coloc := !is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0.5]
atlas[, is_seed := is_deg & has_coloc]
cat("DEGs:", sum(atlas$is_deg, na.rm = TRUE), "\n")
cat("COLOC PP4>0.5 (SuSiE):", sum(atlas$has_coloc, na.rm = TRUE), "\n")
cat("Seeds (DEG + COLOC):", sum(atlas$is_seed, na.rm = TRUE), "\n")

# Standardize column names: atlas uses ensembl_id / human_symbol
# Rename for consistency throughout script
setnames(atlas, "ensembl_id", "gene", skip_absent = TRUE)
setnames(atlas, "human_symbol", "symbol", skip_absent = TRUE)

gene_info <- atlas[, .(gene, symbol, gene_biotype)]

# ---------------------------------------------------------------------------
# 2. Load expression data for coexpression layer
# ---------------------------------------------------------------------------
cat("\n--- Loading expression data ---\n")
dge_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
dge <- readRDS(dge_path)
cat("DGEList:", nrow(dge), "genes x", ncol(dge), "samples\n")

# Compute logCPM
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# Restrict to top 15K genes by mean expression to keep network tractable
gene_means <- rowMeans(logcpm)
keep_expr <- names(sort(gene_means, decreasing = TRUE))[1:min(15000, length(gene_means))]
logcpm_sub <- logcpm[keep_expr, ]
cat("Expression subset:", nrow(logcpm_sub), "genes\n")

# Gene IDs in expression matrix — strip version
expr_genes <- sub("\\.\\d+$", "", rownames(logcpm_sub))
names(expr_genes) <- rownames(logcpm_sub)

# ---------------------------------------------------------------------------
# 3. Build Layer 1: STRING PPI (score > 700)
# ---------------------------------------------------------------------------
cat("\n--- Building PPI layer ---\n")
ppi_file <- file.path(BASE, "data/string_ppi/9606.protein.links.v12.0.txt.gz")
info_file <- file.path(BASE, "data/string_ppi/9606.protein.info.v12.0.txt.gz")

ppi <- fread(ppi_file)
ppi_info <- fread(info_file)
cat("Raw STRING edges:", nrow(ppi), "\n")

# Filter high-confidence
ppi <- ppi[combined_score > 700]
cat("Score>700 edges:", nrow(ppi), "\n")

# Map ENSP to gene symbol
ppi_info[, ensp := sub("^9606\\.", "", `#string_protein_id`)]
ensp2sym <- setNames(ppi_info$preferred_name, ppi_info$ensp)

ppi[, p1 := sub("^9606\\.", "", protein1)]
ppi[, p2 := sub("^9606\\.", "", protein2)]
ppi[, sym1 := ensp2sym[p1]]
ppi[, sym2 := ensp2sym[p2]]

# Drop unmapped
ppi <- ppi[!is.na(sym1) & !is.na(sym2)]
cat("PPI edges with symbols:", nrow(ppi), "\n")

# Map symbols to ENSG IDs from our atlas
sym2ensg <- setNames(atlas$gene, atlas$symbol)
# Strip version from atlas gene IDs
sym2ensg_clean <- sub("\\.\\d+$", "", sym2ensg)
names(sym2ensg_clean) <- names(sym2ensg)

ppi[, g1 := sym2ensg_clean[sym1]]
ppi[, g2 := sym2ensg_clean[sym2]]
ppi <- ppi[!is.na(g1) & !is.na(g2) & g1 != g2]
cat("PPI edges mapped to ENSG:", nrow(ppi), "\n")

# Normalize score to [0,1]
ppi[, w := combined_score / 1000]
ppi_edges <- ppi[, .(from = g1, to = g2, weight = w)]

# ---------------------------------------------------------------------------
# 4. Build Layer 2: Coexpression (top 0.1% correlations)
# ---------------------------------------------------------------------------
cat("\n--- Building coexpression layer ---\n")

# Use expr_genes for mapping
# We'll compute correlations in blocks to manage memory
n_genes <- nrow(logcpm_sub)
cat("Computing correlation matrix for", n_genes, "genes...\n")

# Compute full correlation matrix (15K x 15K manageable at ~1.8GB)
cor_mat <- cor(t(logcpm_sub), method = "pearson")

# Set diagonal and lower triangle to NA to avoid duplicates
cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA

# Find top 0.1% threshold
all_cors <- as.vector(cor_mat)
all_cors <- all_cors[!is.na(all_cors)]
thresh <- quantile(abs(all_cors), probs = 0.999)
cat("Coexpression threshold (top 0.1%):", round(thresh, 3), "\n")

# Extract edges above threshold
idx <- which(abs(cor_mat) >= thresh, arr.ind = TRUE)
cat("Coexpression edges above threshold:", nrow(idx), "\n")

# Map to ENSG IDs (versioned rownames -> unversioned)
rn <- rownames(logcpm_sub)
coexpr_edges <- data.table(
  from = sub("\\.\\d+$", "", rn[idx[, 1]]),
  to   = sub("\\.\\d+$", "", rn[idx[, 2]]),
  weight = abs(cor_mat[idx])
)
coexpr_edges <- coexpr_edges[from != to]
cat("Coexpression edges mapped:", nrow(coexpr_edges), "\n")

# Free memory
rm(cor_mat, all_cors, idx)
gc()

# ---------------------------------------------------------------------------
# 5. Build Layer 3: Regulatory (WGCNA module co-membership)
# ---------------------------------------------------------------------------
cat("\n--- Building regulatory layer ---\n")

# Try to use WGCNA module assignments
wgcna_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier/wgcna_eigengenes.rds")

reg_edges <- data.table()
if (file.exists(wgcna_path)) {
  eigengenes <- readRDS(wgcna_path)
  cat("WGCNA eigengenes loaded:", ncol(eigengenes), "modules\n")

  # We need the module assignments — check if a net object or color vector exists
  wgcna_net_path <- sub("eigengenes", "net", wgcna_path)
  wgcna_colors_path <- sub("eigengenes", "module_colors", wgcna_path)

  # Try to reconstruct from TOM or from the expression data
  # Since we have eigengenes but not module assignments directly,

  # compute module membership via correlation with eigengenes
  # and assign each gene to its best-correlated module

  # Subset to genes in eigengene matrix samples
  common_samples <- intersect(colnames(logcpm_sub), rownames(eigengenes))
  if (length(common_samples) > 50) {
    cat("Computing module memberships from eigengene correlations...\n")
    cat("Common samples:", length(common_samples), "\n")

    me <- as.matrix(eigengenes[common_samples, , drop = FALSE])
    expr_for_mm <- t(logcpm_sub[, common_samples])

    # Correlate each gene with each module eigengene
    mm <- cor(expr_for_mm, me, method = "pearson")

    # Assign each gene to module with highest absolute correlation
    best_module <- apply(abs(mm), 1, which.max)
    best_cor <- apply(abs(mm), 1, max)

    # Only keep genes with |cor| > 0.3 to a module
    assigned <- data.table(
      gene = sub("\\.\\d+$", "", names(best_module)),
      module = best_module,
      mm_cor = best_cor
    )
    assigned <- assigned[mm_cor > 0.3]
    cat("Genes assigned to modules (|cor|>0.3):", nrow(assigned), "\n")

    # Create edges between genes in the same module (weighted by product of memberships)
    # To keep tractable, limit to top 200 genes per module
    modules <- unique(assigned$module)
    edge_list <- list()

    for (m in modules) {
      mg <- assigned[module == m]
      if (nrow(mg) > 200) mg <- mg[order(-mm_cor)][1:200]
      if (nrow(mg) < 2) next

      # Create pairwise edges
      pairs <- combn(seq_len(nrow(mg)), 2)
      edge_list[[length(edge_list) + 1]] <- data.table(
        from = mg$gene[pairs[1, ]],
        to   = mg$gene[pairs[2, ]],
        weight = mg$mm_cor[pairs[1, ]] * mg$mm_cor[pairs[2, ]]
      )
    }

    reg_edges <- rbindlist(edge_list)
    reg_edges <- reg_edges[from != to]
    cat("Regulatory edges (WGCNA co-module):", nrow(reg_edges), "\n")

    rm(me, expr_for_mm, mm, assigned, edge_list)
    gc()
  } else {
    cat("WARNING: Insufficient sample overlap with eigengenes, skipping L3\n")
  }
} else {
  cat("WARNING: WGCNA eigengenes not found, skipping regulatory layer\n")
}

# ---------------------------------------------------------------------------
# 6. Build unified gene universe and multiplex adjacency
# ---------------------------------------------------------------------------
cat("\n--- Building multiplex network ---\n")

# Collect all genes across layers
all_genes <- unique(c(
  ppi_edges$from, ppi_edges$to,
  coexpr_edges$from, coexpr_edges$to,
  reg_edges$from, reg_edges$to
))
cat("Total genes in network:", length(all_genes), "\n")

# Also add seed genes that might not be in any edge list
seed_genes_ensg <- sub("\\.\\d+$", "", atlas[is_seed == TRUE]$gene)
missing_seeds <- setdiff(seed_genes_ensg, all_genes)
cat("Seed genes missing from network:", length(missing_seeds), "(will be added as isolated)\n")
all_genes <- unique(c(all_genes, seed_genes_ensg))

n <- length(all_genes)
gene_idx <- setNames(seq_along(all_genes), all_genes)
cat("Final gene universe:", n, "genes\n")

# Function to build sparse adjacency from edge list
build_adj <- function(edges, gene_idx, n) {
  if (nrow(edges) == 0) return(sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n)))

  i_idx <- gene_idx[edges$from]
  j_idx <- gene_idx[edges$to]
  valid <- !is.na(i_idx) & !is.na(j_idx)

  # Symmetric
  A <- sparseMatrix(
    i = c(i_idx[valid], j_idx[valid]),
    j = c(j_idx[valid], i_idx[valid]),
    x = c(edges$weight[valid], edges$weight[valid]),
    dims = c(n, n)
  )
  return(A)
}

A_ppi    <- build_adj(ppi_edges, gene_idx, n)
A_coexpr <- build_adj(coexpr_edges, gene_idx, n)
A_reg    <- build_adj(reg_edges, gene_idx, n)

cat("PPI adjacency: nnz =", nnzero(A_ppi), "\n")
cat("Coexpression adjacency: nnz =", nnzero(A_coexpr), "\n")
cat("Regulatory adjacency: nnz =", nnzero(A_reg), "\n")

# ---------------------------------------------------------------------------
# 7. Column-normalize each layer and combine (supra-adjacency)
# ---------------------------------------------------------------------------
cat("\n--- Normalizing and combining layers ---\n")

col_normalize <- function(A) {
  cs <- colSums(A)
  cs[cs == 0] <- 1  # avoid division by zero
  # Divide each column by its sum
  A_norm <- A %*% Diagonal(x = 1 / cs)
  return(A_norm)
}

# Normalize each layer
M_ppi    <- col_normalize(A_ppi)
M_coexpr <- col_normalize(A_coexpr)
M_reg    <- col_normalize(A_reg)

# Count active layers per gene (for weighting)
n_layers <- 3
# Simple average across layers (equal weighting)
# Genes present in multiple layers get information from all
M <- (M_ppi + M_coexpr + M_reg) / n_layers

cat("Combined transition matrix: dim =", dim(M), "\n")
cat("Nonzeros in M:", nnzero(M), "\n")

# Clean up individual matrices
rm(A_ppi, A_coexpr, A_reg, M_ppi, M_coexpr, M_reg)
rm(ppi_edges, coexpr_edges, reg_edges)
gc()

# ---------------------------------------------------------------------------
# 8. Random Walk with Restart
# ---------------------------------------------------------------------------
cat("\n--- Running RWR ---\n")

r <- 0.7     # restart probability
tol <- 1e-8  # convergence tolerance
max_iter <- 500

# Initial probability vector: uniform over seeds
seed_idx <- gene_idx[intersect(seed_genes_ensg, names(gene_idx))]
cat("Seeds in network:", length(seed_idx), "of", sum(atlas$is_seed, na.rm = TRUE), "\n")

p0 <- rep(0, n)
p0[seed_idx] <- 1 / length(seed_idx)

p <- p0
for (iter in seq_len(max_iter)) {
  p_new <- (1 - r) * as.vector(M %*% p) + r * p0

  # Normalize to sum to 1
  p_new <- p_new / sum(p_new)

  delta <- max(abs(p_new - p))
  if (iter %% 50 == 0 || delta < tol) {
    cat(sprintf("  Iter %d: max delta = %.2e\n", iter, delta))
  }
  p <- p_new

  if (delta < tol) {
    cat(sprintf("  Converged at iteration %d (delta = %.2e)\n", iter, delta))
    break
  }
}

if (iter == max_iter) cat("  WARNING: Did not converge in", max_iter, "iterations\n")

# ---------------------------------------------------------------------------
# 9. Assemble results
# ---------------------------------------------------------------------------
cat("\n--- Assembling results ---\n")

rwr_dt <- data.table(
  gene = all_genes,
  rwr_score = p
)

# Merge with atlas for annotations
# Atlas gene IDs are already unversioned (after rename above)
atlas_slim <- atlas[, .(gene,
                        symbol, gene_biotype,
                        dream_logFC, dream_padj,
                        coloc_susie_best_pp4,
                        is_deg, has_coloc, is_seed,
                        is_conserved)]

# Handle potential duplicates
atlas_slim <- atlas_slim[!duplicated(gene)]

rwr_dt <- merge(rwr_dt, atlas_slim, by = "gene", all.x = TRUE)
rwr_dt[, rwr_rank := frank(-rwr_score, ties.method = "min")]

# Biotype summary
cat("\nRWR score summary by biotype:\n")
rwr_dt[!is.na(gene_biotype), .(
  n = .N,
  mean_score = mean(rwr_score),
  median_score = median(rwr_score),
  top100 = sum(rwr_rank <= 100)
), by = gene_biotype][order(-mean_score)][1:10] |> print()

# Top lncRNAs
cat("\nTop 20 lncRNAs by RWR score:\n")
rwr_dt[gene_biotype == "lncRNA"][order(rwr_rank)][1:20,
  .(rwr_rank, symbol, rwr_score, is_deg, dream_padj, coloc_susie_best_pp4)] |> print()

# Save
fwrite(rwr_dt[order(rwr_rank)], out_csv)
cat("\nSaved:", out_csv, "\n")

# ---------------------------------------------------------------------------
# 10. Enrichment tests
# ---------------------------------------------------------------------------
cat("\n--- Enrichment tests ---\n")

# Top 500 RWR genes
top500 <- rwr_dt[rwr_rank <= 500]$gene

# Background = genes in the network (not full atlas, which includes genes absent from network)
bg_genes <- rwr_dt$gene

# DEG enrichment (restricted to network genes)
deg_genes <- intersect(rwr_dt[is_deg == TRUE]$gene, bg_genes)
a <- length(intersect(top500, deg_genes))
b <- length(setdiff(top500, deg_genes))
c_val <- length(setdiff(deg_genes, top500))
d <- length(setdiff(bg_genes, union(top500, deg_genes)))
ft_deg <- fisher.test(matrix(c(a, b, c_val, d), nrow = 2), alternative = "greater")
cat(sprintf("DEG enrichment: OR=%.2f, p=%.2e (%d/%d in top500)\n",
            ft_deg$estimate, ft_deg$p.value, a, a + b))

# COLOC enrichment (restricted to network genes)
coloc_genes <- intersect(rwr_dt[has_coloc == TRUE]$gene, bg_genes)
a2 <- length(intersect(top500, coloc_genes))
b2 <- length(setdiff(top500, coloc_genes))
c2 <- length(setdiff(coloc_genes, top500))
d2 <- length(setdiff(bg_genes, union(top500, coloc_genes)))
ft_coloc <- fisher.test(matrix(c(a2, b2, c2, d2), nrow = 2), alternative = "greater")
cat(sprintf("COLOC enrichment: OR=%.2f, p=%.2e (%d/%d in top500)\n",
            ft_coloc$estimate, ft_coloc$p.value, a2, a2 + b2))

# Conserved enrichment (restricted to network genes)
cc_genes <- intersect(rwr_dt[is_conserved == TRUE]$gene, bg_genes)
cc_genes <- cc_genes[!is.na(cc_genes)]
a3 <- length(intersect(top500, cc_genes))
b3 <- length(setdiff(top500, cc_genes))
c3 <- length(setdiff(cc_genes, top500))
d3 <- length(setdiff(bg_genes, union(top500, cc_genes)))
ft_cc <- fisher.test(matrix(c(a3, b3, c3, d3), nrow = 2), alternative = "greater")
cat(sprintf("Conserved enrichment: OR=%.2f, p=%.2e (%d/%d in top500)\n",
            ft_cc$estimate, ft_cc$p.value, a3, a3 + b3))

# lncRNA DEG enrichment among RWR-scored lncRNAs
lncrna_in_net <- rwr_dt[gene_biotype == "lncRNA" & !is.na(rwr_score)]
lnc_top <- lncrna_in_net[rwr_rank <= quantile(rwr_rank, 0.25)]$gene
lnc_deg <- lncrna_in_net[is_deg == TRUE]$gene
a4 <- length(intersect(lnc_top, lnc_deg))
b4 <- length(setdiff(lnc_top, lnc_deg))
c4 <- length(setdiff(lnc_deg, lnc_top))
d4 <- nrow(lncrna_in_net) - length(union(lnc_top, lnc_deg))
if (d4 > 0) {
  ft_lnc <- fisher.test(matrix(c(a4, b4, c4, d4), nrow = 2), alternative = "greater")
  cat(sprintf("lncRNA top25%% vs DEG: OR=%.2f, p=%.2e (%d/%d)\n",
              ft_lnc$estimate, ft_lnc$p.value, a4, a4 + b4))
}

enrichment_results <- data.table(
  test = c("DEG", "COLOC", "Conserved"),
  overlap = c(a, a2, a3),
  top500 = 500,
  OR = c(ft_deg$estimate, ft_coloc$estimate, ft_cc$estimate),
  pvalue = c(ft_deg$p.value, ft_coloc$p.value, ft_cc$p.value)
)
cat("\nEnrichment summary:\n")
print(enrichment_results)

# ---------------------------------------------------------------------------
# 11. Figures (3 panels)
# ---------------------------------------------------------------------------
cat("\n--- Generating figures ---\n")

# Biotype mapping for clean labels
rwr_dt[, biotype_label := fcase(
  gene_biotype == "protein_coding", "Protein-coding",
  gene_biotype == "lncRNA", "lncRNA",
  gene_biotype %in% c("miRNA", "snoRNA", "snRNA"), "Small ncRNA",
  !is.na(gene_biotype), "Other",
  default = "Unknown"
)]

# Panel (a): RWR score distribution by biotype
p_dist <- ggplot(rwr_dt[biotype_label %in% c("Protein-coding", "lncRNA", "Small ncRNA", "Other")],
       aes(x = biotype_label, y = log10(rwr_score + 1e-10), fill = biotype_label)) +
  geom_violin(scale = "width", alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3, alpha = 0.7) +
  scale_fill_manual(values = c(
    "Protein-coding" = masld_colors$deg,
    "lncRNA" = masld_colors$up,
    "Small ncRNA" = masld_colors$twas,
    "Other" = masld_colors$ns
  )) +
  labs(x = NULL, y = expression(log[10](RWR~score)),
       title = "RWR score by gene biotype") +
  theme_masld() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

# Panel (b): Top 20 lncRNAs by RWR score
top20_lnc <- rwr_dt[gene_biotype == "lncRNA"][order(rwr_rank)][1:20]
top20_lnc[, deg_label := fifelse(is_deg == TRUE, "DEG", "Not DEG")]
top20_lnc[, symbol_display := fifelse(is.na(symbol) | symbol == "", gene, symbol)]

p_lnc <- ggplot(top20_lnc, aes(x = reorder(symbol_display, rwr_score), y = rwr_score, fill = deg_label)) +
  geom_col(alpha = 0.9) +
  coord_flip() +
  scale_fill_manual(values = c("DEG" = masld_colors$up, "Not DEG" = masld_colors$ns),
                    name = NULL) +
  labs(x = NULL, y = "RWR score",
       title = "Top 20 lncRNAs by RWR score") +
  theme_masld() +
  theme(legend.position = c(0.7, 0.3))

# Panel (c): Enrichment barplot
enrichment_results[, neg_log10p := -log10(pvalue)]
enrichment_results[, sig_label := fifelse(pvalue < 0.05, "*", "")]
enrichment_results[pvalue < 0.01, sig_label := "**"]
enrichment_results[pvalue < 0.001, sig_label := "***"]

p_enrich <- ggplot(enrichment_results, aes(x = test, y = neg_log10p, fill = test)) +
  geom_col(alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text(aes(label = sprintf("OR=%.1f\n%s", OR, sig_label)),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(
    "DEG" = masld_colors$deg,
    "COLOC" = masld_colors$gwas,
    "Conserved" = masld_colors$conserved
  )) +
  labs(x = NULL, y = expression(-log[10](p)),
       title = "RWR top-500 enrichment") +
  theme_masld() +
  theme(legend.position = "none")

# Combine panels
pdf(out_pdf, width = 12, height = 4.5)
gridExtra::grid.arrange(p_dist, p_lnc, p_enrich, ncol = 3,
                        widths = c(1.1, 1.2, 1))
dev.off()
cat("Saved figure:", out_pdf, "\n")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Genes scored:", nrow(rwr_dt), "\n")
cat("lncRNAs in network:", sum(rwr_dt$gene_biotype == "lncRNA", na.rm = TRUE), "\n")
cat("lncRNAs in top 500:", sum(rwr_dt$gene_biotype == "lncRNA" & rwr_dt$rwr_rank <= 500, na.rm = TRUE), "\n")
cat("Seeds used:", length(seed_idx), "\n")
cat("RWR iterations:", iter, "\n")
cat("Network layers: PPI + Coexpression + Regulatory (WGCNA)\n")
cat("Output CSV:", out_csv, "\n")
cat("Output PDF:", out_pdf, "\n")
cat("Done:", format(Sys.time()), "\n")
