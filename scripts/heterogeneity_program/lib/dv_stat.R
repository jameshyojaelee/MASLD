# ─────────────────────────────────────────────────────────────────────────────
# T1 / Phase-0c — limma differential-variability statistic + fold-internal callable
#
# `diffVar` IS limma on absolute deviations (missMethyl::varFit logic), so we
# stay in the limma family with NO new engine. Brown-Forsythe form: centre each
# gene on its per (block×group) MEDIAN (robust), then run limma `lmFit`+`eBayes`
# on the |deviations|; the group coefficient's moderated t is the DV statistic.
# Counts → voom log-CPM with `dataset` in the design (NEVER removeBatchEffect —
# GATE-N8). Quality weights optional (user: stay in limma family, not bound to
# voom-QW). Positive DV t = MORE variable in the non-reference group (disease).
#
# `dv_screen(train_idx)` is the GATE-A fold-internal callable handed to T2/T3:
# stateless, recomputes filterByExpr+voom+DV from scratch on the given samples.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(limma); library(edgeR) })

.have_ms <- requireNamespace("matrixStats", quietly = TRUE)
.rowMedians <- function(m) {
  if (ncol(m) == 1L) return(as.numeric(m))
  if (.have_ms) matrixStats::rowMedians(m) else apply(m, 1L, median)
}

# Brown-Forsythe / diffVar moderated DV statistic.
#  E       : G x N voom log-CPM
#  weights : G x N voom precision weights (or NULL)
#  group   : length-N factor; level 1 = reference (control)
#  block   : length-N blocking factor (dataset)
# returns list(t, logfc, p) for the group coefficient on |deviations|.
dv_diffvar <- function(E, weights, group, block) {
  group <- droplevels(as.factor(group)); block <- droplevels(as.factor(block))
  cell  <- droplevels(interaction(block, group))
  z <- E
  for (cl in levels(cell)) {
    idx <- which(cell == cl)
    med <- .rowMedians(E[, idx, drop = FALSE])
    z[, idx] <- abs(E[, idx, drop = FALSE] - med)   # absolute deviation from cell median
  }
  design <- model.matrix(~ block + group)
  gcol <- tail(grep("^group", colnames(design)), 1L)   # last group column = DV contrast
  fit <- lmFit(z, design, weights = weights)
  fit <- eBayes(fit, trend = TRUE)                      # trend handles |dev| mean-variance
  list(t = fit$t[, gcol], logfc = fit$coefficients[, gcol], p = fit$p.value[, gcol])
}

# Fold-internal DV screen (GATE-A). dge = DGEList (raw counts); group/block length = ncol(dge).
# gene_set: optional fixed gene vector (for permutation — keep the universe constant
# across replicates while refitting voom+DV).  qw: include quality weights.
dv_screen <- function(dge, group, block, train_idx = seq_len(ncol(dge)),
                      gene_set = NULL, qw = FALSE) {
  d <- dge[, train_idx, keep.lib.sizes = FALSE]
  g <- droplevels(as.factor(group[train_idx])); b <- droplevels(as.factor(block[train_idx]))
  if (is.null(gene_set)) {
    keep <- filterByExpr(d, group = g)
    d <- d[keep, , keep.lib.sizes = FALSE]
  } else {
    d <- d[rownames(d) %in% gene_set, , keep.lib.sizes = FALSE]
  }
  des <- model.matrix(~ b + g)
  v <- if (qw) voomWithQualityWeights(d, des) else voom(d, des)
  res <- dv_diffvar(v$E, v$weights, g, b)
  data.frame(gene = rownames(d), dv_t = res$t, dv_logfc = res$logfc,
             dv_p = res$p, AveExpr = rowMeans(v$E), stringsAsFactors = FALSE)
}

# Count-native NB confirmer (edgeR): per-group tagwise dispersion log-ratio.
# Concordance arm only (NOT in the permutation significance). Positive = more
# overdispersed in the non-reference (disease) group.
dv_nb_confirmer <- function(dge, group, block, gene_set) {
  g <- droplevels(as.factor(group))
  d <- dge[rownames(dge) %in% gene_set, , keep.lib.sizes = FALSE]
  lvl <- levels(g)
  disp <- lapply(lvl, function(L) {
    dl <- d[, g == L, keep.lib.sizes = FALSE]
    dl <- estimateDisp(dl, design = model.matrix(~ droplevels(block[g == L])), robust = TRUE)
    dl$tagwise.dispersion
  })
  data.frame(gene = rownames(d),
             nb_logratio = log(pmax(disp[[2]], 1e-8) / pmax(disp[[1]], 1e-8)),
             stringsAsFactors = FALSE)
}
