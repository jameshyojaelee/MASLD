#!/usr/bin/env Rscript
# 03: the co-expression graph the modules and their competitors both come from.
#
# Built on the discovery cohorts, on the shared universe, with no outcome read.
# Two choices matter and are made here rather than left implicit.
#
# Genes are z-scored WITHIN cohort before the correlation. Pooling five cohorts
# raw would let a between-cohort mean shift create correlation between any two
# genes that happen to differ between cohorts, which is a batch structure, not
# liver co-expression.
#
# Only POSITIVE correlations become edges. A module here is scored as a mean of
# its members, so a member anti-correlated with the rest cancels signal rather
# than adding it; an edge that cannot contribute to the score should not shape
# the partition.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))
suppressPackageStartupMessages({
  library(edgeR)
  library(igraph)
})

contract <- cam_contract()
out <- cam_dir("graph")
cam_assert(file.exists(file.path(cam_out_root(), "universe", "READY")),
           "Run 02_build_universe.R first")

universe <- fread(file.path(cam_out_root(), "universe", "universe.tsv"))$gene_symbol
manifest <- fread(cam_input("bulk_manifest", contract))
cam_assert(nrow(manifest) == 844L, paste0("Expected 844 manifest rows, got ", nrow(manifest)))

cam_say("loading five-cohort DGEList")
dge <- readRDS(cam_input("bulk_dge", contract))
cam_assert(all(manifest$sample_id %in% colnames(dge)),
           "Manifest samples are not all present in the DGEList")
dge <- dge[, manifest$sample_id]
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# ENSG to symbol from the substrate's own results file, so the module gene
# space and the differential-expression gene space cannot diverge.
map <- unique(fread(cam_input("bulk_stage_results", contract),
                    select = c("gene_id_versioned", "gene_name")))
sym <- map[match(rownames(logcpm), gene_id_versioned), gene_name]
keep <- !is.na(sym) & sym %in% universe
logcpm <- logcpm[keep, , drop = FALSE]
sym <- sym[keep]

# Several ENSG ids can carry one symbol. Averaging on the log scale keeps the
# collapsed value on the scale the z-score expects.
cam_say("collapsing ", nrow(logcpm), " ENSG rows to symbols")
expr <- rowsum(logcpm, group = sym, reorder = TRUE) /
  as.vector(table(sym)[sort(unique(sym))])
expr <- expr[intersect(universe, rownames(expr)), , drop = FALSE]
cam_say("expression matrix: ", nrow(expr), " symbols x ", ncol(expr), " participants")
cam_assert(nrow(expr) >= contract$universe$min_universe_size,
           "Universe shrank below the floor during symbol collapse")

cam_say("z-scoring within cohort")
z <- expr
for (ds in unique(manifest$dataset)) {
  j <- which(manifest$dataset == ds)
  block <- expr[, j, drop = FALSE]
  mu <- rowMeans(block)
  sdv <- apply(block, 1L, stats::sd)
  sdv[!is.finite(sdv) | sdv == 0] <- NA_real_
  z[, j] <- (block - mu) / sdv
}
# A gene with no variance in one cohort carries no ordering there; dropping it
# from the graph is more honest than imputing a zero row.
ok <- stats::complete.cases(z)
cam_say("genes with finite within-cohort z in every cohort: ", sum(ok), " of ", nrow(z))
z <- z[ok, , drop = FALSE]
cam_write_rds(z, file.path(out, "cohort_z_expression.rds"))
# Mean log-CPM per gene BEFORE centring. The expression-matched null bins on
# this; binning the within-cohort z (v1) binned numerical zeros.
cam_write_tsv(data.table(gene_symbol = rownames(expr), mean_logcpm = rowMeans(expr)),
              file.path(out, "gene_mean_logcpm.tsv"))

cam_say("Spearman correlation across ", ncol(z), " participants")
rho <- stats::cor(t(z), method = "spearman")
diag(rho) <- 0
# Unclipped correlations are kept so the coherence of a module can be compared
# with the coherence of its size-matched competitors on the same substrate.
cam_write_rds(rho, file.path(out, "spearman_rho_unclipped.rds"))
rho[rho < 0] <- 0

# Neighbourhood size. The prespecified k = 10 left 38 percent of the universe
# with no edge at all on a pool this size, so k is selected from a grid by a
# connectivity rule that reads no outcome. See GRAPH_K_AMENDMENT.md.
k_grid <- contract$graph$k_grid
build_mutual <- function(k) {
  nn <- t(apply(rho, 1L, function(r) order(r, decreasing = TRUE)[seq_len(k)]))
  flag <- matrix(FALSE, nrow(rho), nrow(rho))
  for (i in seq_len(nrow(rho))) flag[i, nn[i, ]] <- TRUE
  m <- flag & t(flag)
  m[rho <= 0] <- FALSE
  m
}
scan_rows <- list()
graphs <- list()
for (k in k_grid) {
  m <- build_mutual(k)
  gk <- igraph::graph_from_adjacency_matrix(m, mode = "undirected", diag = FALSE)
  igraph::V(gk)$name <- rownames(z)
  ck <- igraph::components(gk)
  deg <- igraph::degree(gk)
  scan_rows[[length(scan_rows) + 1L]] <- data.table(
    k = k, n_edges = igraph::gsize(gk), n_components = ck$no,
    largest_component_size = max(ck$csize),
    largest_component_fraction = max(ck$csize) / nrow(rho),
    median_degree = stats::median(deg), n_isolated = sum(deg == 0),
    isolated_fraction = mean(deg == 0)
  )
  graphs[[as.character(k)]] <- gk
  cam_say("k=", k, " edges ", igraph::gsize(gk),
          " largest ", round(max(ck$csize) / nrow(rho), 3),
          " median degree ", stats::median(deg),
          " isolated ", sum(deg == 0))
}
scan <- rbindlist(scan_rows)
cam_write_tsv(scan, file.path(out, "k_selection_scan.tsv"))

eligible <- scan[largest_component_fraction >= contract$graph$k_min_main_component_fraction_for_selection &
                   median_degree >= contract$graph$k_min_median_degree]
k_selected <- if (nrow(eligible)) min(eligible$k) else max(k_grid)
sparse_fallback <- nrow(eligible) == 0L
cam_say("selected k = ", k_selected, if (sparse_fallback) " (fallback: no k met the rule)" else "")

g <- graphs[[as.character(k_selected)]]
mutual <- build_mutual(k_selected)
igraph::E(g)$weight <- rho[igraph::as_edgelist(g, names = FALSE)]
comp <- igraph::components(g)
n <- nrow(rho)
main_fraction <- max(comp$csize) / n
cam_say("nodes ", n, " edges ", igraph::gsize(g),
        " components ", comp$no, " largest fraction ", round(main_fraction, 4))
cam_assert(main_fraction >= contract$graph$min_main_component_fraction,
  paste0("Largest component holds ", round(main_fraction, 3),
         " of the universe, below the prespecified floor of ",
         contract$graph$min_main_component_fraction, "; the analysis stops here"))

cam_write_rds(g, file.path(out, "coexpression_graph.rds"))
adjacency <- igraph::as_adj_list(g, mode = "all")
adjacency <- lapply(adjacency, as.integer)
cam_write_rds(adjacency, file.path(out, "adjacency_list.rds"))
cam_write_tsv(
  data.table(gene_symbol = igraph::V(g)$name, degree = igraph::degree(g),
             component = comp$membership),
  file.path(out, "graph_nodes.tsv")
)
cam_write_json(list(
  n_nodes = n, n_edges = igraph::gsize(g), n_components = comp$no,
  largest_component_size = max(comp$csize),
  largest_component_fraction = main_fraction,
  k_selected = k_selected, k_grid = k_grid,
  k_selection_was_fallback = sparse_fallback,
  edge_rule = contract$graph$edge_rule,
  median_degree = stats::median(igraph::degree(g)),
  n_isolated = sum(igraph::degree(g) == 0)
), file.path(out, "graph_summary.json"))

writeLines("graph built", file.path(out, "READY"))
cam_say("03 complete")
