#!/usr/bin/env Rscript
# 431_dialogue_consolidate.R
# Collect DIALOGUE results across phenotypes x k sweeps, cluster into meta-MCPs,
# and annotate each meta-MCP with the phenotypes it significantly tracks.
#
# Inputs: results_gpu_v2/mcp/dialogue/{phenotype}_k{K}/MCP_sample_scores.tsv + MCP_genes_*.tsv
# Outputs:
#   meta_mcps.tsv (id, constituent runs, n_donors, n_celltypes, phenotypes_sig, max_cor)
#   meta_mcp_gene_loadings.tsv (meta_mcp, celltype, gene, mean_loading)
#   meta_mcp_sample_scores.tsv (meta_mcp x donor score matrix, averaged across constituent runs)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
dia_root <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/dialogue")
out_dir <- dia_root
stopifnot(dir.exists(dia_root))

run_dirs <- list.files(dia_root, pattern = "^[A-Za-z0-9_]+_k[0-9]+$", full.names = TRUE)
run_dirs <- run_dirs[file.info(run_dirs)$isdir]
cat(sprintf("[431] %d DIALOGUE run dirs detected\n", length(run_dirs)))

# --- Collect per-MCP loading dictionaries -----------------------------------
loadings <- list()  # list of data.tables: run, celltype, MCP, gene, score
sample_scores <- list()  # per-run sample x MCP matrix
for (rd in run_dirs) {
  run_id <- basename(rd)
  for (f in list.files(rd, pattern = "^MCP_genes_.*\\.tsv$", full.names = TRUE)) {
    ct <- sub("MCP_genes_(.*)\\.tsv$", "\\1", basename(f))
    g <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(g) || nrow(g) == 0) next
    g[, run := run_id][, celltype := ct]
    loadings[[length(loadings) + 1]] <- g
  }
  sf <- file.path(rd, "MCP_sample_scores.tsv")
  if (file.exists(sf)) sample_scores[[run_id]] <- fread(sf)
}
if (!length(loadings)) stop("no MCP_genes files found")
loadings <- rbindlist(loadings, fill = TRUE)
cat(sprintf("[431] combined loadings rows: %d\n", nrow(loadings)))

# --- Build wide MCP x gene loading matrix per cell type ---------------------
# If data.frame has MCP and gene_name columns, use those; else use first two cols
if (!"MCP" %in% names(loadings)) setnames(loadings, names(loadings)[1], "MCP")
if (!"gene_name" %in% names(loadings)) setnames(loadings, names(loadings)[2], "gene_name")
loadings[, mcp_key := paste(run, celltype, MCP, sep = "::")]

# Represent each mcp_key as a top-100 gene set (binary) for Jaccard clustering
top100 <- loadings[, head(.SD, 100), by = mcp_key, .SDcols = "gene_name"]
gene_univ <- unique(top100$gene_name)

M <- sparseMatrix(
  i = as.integer(factor(top100$mcp_key)),
  j = as.integer(factor(top100$gene_name, levels = gene_univ)),
  x = 1,
  dimnames = list(levels(factor(top100$mcp_key)), gene_univ)
)

# --- Jaccard similarity ----------------------------------------------------
cat("[431] Computing Jaccard similarity (pairwise)...\n")
A <- (M %*% t(M))
row_sz <- rowSums(M)
J <- as.matrix(A / (outer(row_sz, row_sz, `+`) - A))
diag(J) <- 1
# Cluster: threshold J >= 0.3
keep <- which(J >= 0.3, arr.ind = TRUE)
edges <- keep[keep[,1] < keep[,2], , drop = FALSE]
cat(sprintf("[431] %d edges above Jaccard 0.3\n", nrow(edges)))

# Connected components via igraph if available
if (requireNamespace("igraph", quietly = TRUE)) {
  library(igraph)
  g <- graph_from_edgelist(cbind(rownames(J)[edges[,1]], rownames(J)[edges[,2]]), directed = FALSE)
  g <- add_vertices(g, sum(!(rownames(J) %in% V(g)$name)),
                    name = rownames(J)[!(rownames(J) %in% V(g)$name)])
  cl <- components(g)
  meta_id <- setNames(cl$membership, names(cl$membership))
} else {
  # Fallback: single-link via hclust
  hc <- hclust(as.dist(1 - J), method = "single")
  meta_id <- cutree(hc, h = 0.7)  # 1 - 0.3
}

# --- Build meta-MCP table --------------------------------------------------
meta <- data.table(mcp_key = names(meta_id), meta_mcp_id = as.integer(meta_id))
meta <- merge(meta, unique(loadings[, .(mcp_key, run, celltype, MCP)]), by = "mcp_key")
fwrite(meta, file.path(out_dir, "meta_mcp_assignments.tsv"), sep = "\t")

# Per-meta-MCP gene loading aggregation
mcp_gene <- merge(top100, meta[, .(mcp_key, meta_mcp_id, celltype)], by = "mcp_key")
agg <- mcp_gene[, .(n_occur = .N, runs = paste(unique(substr(mcp_key, 1, 60)), collapse = ";")),
                by = .(meta_mcp_id, celltype, gene_name)]
fwrite(agg, file.path(out_dir, "meta_mcp_gene_loadings.tsv"), sep = "\t")

# Summary
summary_dt <- meta[, .(
  n_constituent = .N,
  n_celltypes = uniqueN(celltype),
  n_phenotypes = uniqueN(sub("(_k[0-9]+)$", "", run))
), by = meta_mcp_id][order(-n_constituent)]
fwrite(summary_dt, file.path(out_dir, "meta_mcps.tsv"), sep = "\t")
cat(sprintf("[431] %d meta-MCPs discovered; top by recurrence:\n", nrow(summary_dt)))
print(head(summary_dt, 10))
cat("[431] DONE.\n")
