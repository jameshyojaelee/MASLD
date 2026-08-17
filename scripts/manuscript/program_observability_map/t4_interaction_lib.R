#!/usr/bin/env Rscript

# T4 (SYS-I): residual program interaction structure.
#
# THE OBJECT
# For every pair of the 117 frozen programs, the correlation of their
# donor-level scores AFTER the recorded histology has been removed from both.
# A pair that covaries only because both track fibrosis contributes nothing
# here; what is left is coordination among donors with the SAME recorded
# histology. This is not a co-expression network: no edge is defined by gene
# co-expression, no module is discovered, and no new data is used.
#
# THE THREE WAYS THIS ANALYSIS CAN LIE, AND WHAT IS DONE ABOUT EACH
#
# 1. STAGE STRUCTURE LEFT IN THE RESIDUAL. A linear stage term leaves the
#    non-linear part of the stage response in the residual; every program that
#    responds non-linearly then keeps a shared stage-driven component and the
#    pairs covary through it. The axis map shows 23 fibrosis programs whose
#    shape a linear term misses, and a sibling system in this workstream had to
#    retract a result for exactly this reason. The residual model here is
#    therefore SATURATED in fibrosis (5-level factor) and coarsened-saturated in
#    NAS (the not-NASH / borderline / NASH cut, whose raw 0-8 levels are too
#    thin at n = 76). t4_stage_leakage() gates the result on there being no
#    residual stage structure left.
#
# 2. SHARED MEMBERSHIP. Programs share genes: 1,698 of 3,498 mapped symbols sit
#    in more than one program. Two programs sharing genes covary even when
#    nothing biological connects them, because part of each score is literally
#    the same measurement. This is handled twice over:
#      (a) the null permutes donors independently PER GENE, which destroys all
#          gene-gene coordination while preserving each program's weight vector
#          and every shared gene exactly. The null mean of an edge is therefore
#          the correlation shared membership alone produces, and the observed
#          statistic is centered on it;
#      (b) the primary edge statistic is recomputed with the shared genes
#          DELETED FROM BOTH programs. (a) is exact only under gene
#          independence and under-corrects when the shared genes are themselves
#          coordinated; (b) has no such loophole. A pair that cannot give up its
#          shared genes and keep 80 percent of its weight is not testable and is
#          reported as membership-entangled rather than as an edge.
#
# 3. POOLING ACROSS COHORTS. Correlating scores pooled over cohorts manufactures
#    edges out of cohort mean differences. Every correlation here is estimated
#    WITHIN cohort and combined on the Fisher scale.
#
# WHAT MAKES THE RESIDUAL SHORTCUT EXACT
# Residualization and program scoring are both linear in the gene matrix, so the
# studentized residual of any program score equals U %*% w, where U is the
# studentization-scaled residual of the gene z-score matrix and w is the
# program's weight vector. Building U once per cohort makes the leave-shared-
# genes-out recomputation and the permutation null affordable. t4_run_tests.R
# checks U %*% w against hc3_fit_matrix() on the real scoring path.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# The NAS coarsening. Same cut as the rest of the workstream.
t4_nas_bin <- function(x) {
  cut(x, breaks = c(-Inf, 2, 4, Inf), labels = c("le2", "3to4", "ge5"),
      right = TRUE)
}

T4_RESIDUAL_FORMULA <- ~ sex_final + fibrosis_factor + nas_factor

# ---------------------------------------------------------------- construction

# Long membership table -> named list of (symbol, weight) per program, on the
# same definition score_programs() uses: mapped_symbol, original_l1_weight,
# summed over duplicate symbols.
t4_membership_weights <- function(membership, features) {
  m <- copy(membership)
  m[, gene_symbol := toupper(trimws(mapped_symbol))]
  m <- m[!is.na(gene_symbol) & gene_symbol != "" & program_uid %in% features]
  m[, .(weight = sum(as.numeric(original_l1_weight))),
    by = .(feature_id = program_uid, gene_symbol)]
}

# One block per cohort: the gene panel that is variable in this cohort, the
# saturated design, the hat-matrix pieces, and the per-cohort weight matrix.
#
# Weights are renormalized over the genes retained IN THIS COHORT, which is what
# score_programs() does; a program whose retained weight falls below the coverage
# threshold in a cohort is dropped from that cohort rather than scored on a
# fraction of itself.
t4_build_block <- function(symbol_matrix, meta_cohort, weights, features,
                           coverage_threshold) {
  Z <- zscore_rows(symbol_matrix[, meta_cohort$sample_id, drop = FALSE])
  variable <- rownames(Z)[rowSums(is.finite(Z)) == ncol(Z)]
  Z <- Z[variable, , drop = FALSE]

  d <- copy(meta_cohort)
  d[, `:=`(fibrosis_factor = droplevels(factor(fibrosis_stage)),
           nas_factor = droplevels(t4_nas_bin(nas_score)),
           sex_final = droplevels(factor(sex_final)))]
  X <- stats::model.matrix(T4_RESIDUAL_FORMULA, data = d)
  if (qr(X)$rank != ncol(X) || nrow(X) <= ncol(X) + 3L) {
    fail("Saturated residual design not estimable in cohort ", meta_cohort$dataset[[1L]])
  }
  xtx_inv <- solve(crossprod(X))
  leverage <- rowSums((X %*% xtx_inv) * X)

  total <- weights[, .(total = sum(weight)), by = feature_id]
  kept <- weights[gene_symbol %in% variable]
  retained <- kept[, .(retained = sum(weight)), by = feature_id]
  cover <- merge(data.table(feature_id = features), total, by = "feature_id", all.x = TRUE)
  cover <- merge(cover, retained, by = "feature_id", all.x = TRUE)
  cover[is.na(retained), retained := 0]
  cover[, `:=`(cohort_weight_coverage = retained / total,
               cohort_testable = is.finite(retained / total) &
                 retained / total >= coverage_threshold)]

  W <- matrix(0, nrow = length(variable), ncol = length(features),
              dimnames = list(variable, features))
  usable <- cover[cohort_testable == TRUE, feature_id]
  kept <- merge(kept, cover[, .(feature_id, retained)], by = "feature_id")
  kept <- kept[feature_id %in% usable]
  W[cbind(match(kept$gene_symbol, variable), match(kept$feature_id, features))] <-
    kept$weight / kept$retained

  list(cohort = meta_cohort$dataset[[1L]], samples = meta_cohort$sample_id,
       genes = variable, Z = Z, X = X, xtx_inv = xtx_inv,
       dinv = 1 / sqrt(pmax(1 - leverage, 1e-8)),
       n = nrow(X), n_par = ncol(X), df = nrow(X) - ncol(X),
       W = W, coverage = cover,
       meta = d[, .(sample_id, fibrosis_stage, nas_score, sex_final)])
}

# Studentization-scaled residual of the gene z-score matrix. Rows are samples,
# columns genes, so U %*% w is the residual score of any weight vector w.
t4_residual_gene_matrix <- function(block, Z = block$Z) {
  Yt <- t(Z)
  B <- block$xtx_inv %*% crossprod(block$X, Yt)
  (Yt - block$X %*% B) * block$dinv
}

# --------------------------------------------------------------- edge geometry

# Every unordered pair, plus the two membership-overlap measures reported
# alongside every edge. Jaccard is the set measure the contract asks for; the
# weight cosine is the analytic prediction of the correlation shared membership
# induces under gene independence, so it is the one that belongs next to a
# correlation.
t4_pair_table <- function(features, weights) {
  np <- length(features)
  ij <- which(upper.tri(matrix(0, np, np)), arr.ind = TRUE)
  pairs <- data.table(pair_index = seq_len(nrow(ij)),
                      i = ij[, "row"], j = ij[, "col"])
  pairs[, `:=`(program_a = features[i], program_b = features[j])]

  sets <- split(weights$gene_symbol, weights$feature_id)[features]
  wide <- dcast(weights, gene_symbol ~ feature_id, value.var = "weight", fill = 0)
  Wm <- as.matrix(wide[, -1L])
  rownames(Wm) <- wide$gene_symbol
  Wm <- Wm[, features, drop = FALSE]
  Wn <- sweep(Wm, 2L, sqrt(colSums(Wm^2)), "/")
  cosine <- crossprod(Wn)
  l1 <- sweep(Wm, 2L, colSums(Wm), "/")

  pairs[, n_shared := vapply(seq_len(.N), function(k)
    length(intersect(sets[[i[k]]], sets[[j[k]]])), integer(1))]
  pairs[, n_union := vapply(seq_len(.N), function(k)
    length(union(sets[[i[k]]], sets[[j[k]]])), integer(1))]
  pairs[, jaccard := n_shared / n_union]
  pairs[, weight_cosine := cosine[cbind(i, j)]]
  pairs[, shared_l1_weight_a := vapply(seq_len(.N), function(k)
    sum(l1[intersect(sets[[i[k]]], sets[[j[k]]]), i[k]]), numeric(1))]
  pairs[, shared_l1_weight_b := vapply(seq_len(.N), function(k)
    sum(l1[intersect(sets[[i[k]]], sets[[j[k]]]), j[k]]), numeric(1))]
  pairs[]
}

# Sparse leave-shared-genes-out operators, per cohort. For pair (a, b) the
# disjoint score of a is (P_a - U %*% v_a) / (1 - s_a), where v_a is a's weight
# vector masked to the shared genes and s_a is the masked weight. Assembling all
# of them as one sparse matrix turns the whole leave-out arm into two matrix
# products per cohort.
t4_disjoint_operators <- function(block, pairs, features, min_retained_weight) {
  overlapping <- pairs[n_shared > 0L]
  genes <- block$genes
  W <- block$W
  npair <- nrow(overlapping)
  s_a <- numeric(npair); s_b <- numeric(npair)
  rows <- vector("list", npair); cols <- vector("list", npair)
  xa <- vector("list", npair); xb <- vector("list", npair)
  sets <- lapply(seq_along(features), function(k) which(W[, k] != 0))
  ivec <- overlapping$i; jvec <- overlapping$j
  for (k in seq_len(npair)) {
    shared <- intersect(sets[[ivec[k]]], sets[[jvec[k]]])
    if (!length(shared)) next
    rows[[k]] <- shared
    cols[[k]] <- rep.int(k, length(shared))
    xa[[k]] <- W[shared, ivec[k]]
    xb[[k]] <- W[shared, jvec[k]]
    s_a[k] <- sum(xa[[k]]); s_b[k] <- sum(xb[[k]])
  }
  ri <- unlist(rows, use.names = FALSE); ci <- unlist(cols, use.names = FALSE)
  build <- function(xs) {
    x <- unlist(xs, use.names = FALSE)
    if (!length(x)) {
      return(sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                          dims = c(length(genes), npair)))
    }
    sparseMatrix(i = ri, j = ci, x = x, dims = c(length(genes), npair))
  }
  list(pair_index = overlapping$pair_index,
       i = ivec, j = jvec,
       Va = build(xa), Vb = build(xb),
       shared_a = s_a, shared_b = s_b,
       eligible = (1 - s_a) >= min_retained_weight &
         (1 - s_b) >= min_retained_weight)
}

# ------------------------------------------------------------------- one pass

# Everything a single realisation of the data produces: per-cohort Fisher z for
# the full-membership and leave-shared-out edge sets, their inverse-variance
# combination, and the global structure statistics. The null calls this with a
# permuted gene matrix and nothing else changes.
T4_PARTIAL_COMPONENTS <- c(1L, 3L)

t4_one_pass <- function(blocks, pairs, ops, Zlist = NULL, want_matrix = FALSE,
                        n_components = T4_PARTIAL_COMPONENTS) {
  np <- nrow(pairs)
  nc <- length(blocks)
  nf <- ncol(blocks[[1L]]$W)
  zc_full <- matrix(NA_real_, np, nc)
  zc_disj <- matrix(NA_real_, np, nc)
  zc_part <- lapply(n_components, function(k) matrix(NA_real_, np, nc))
  weight_c <- numeric(nc)
  cor_sum <- matrix(0, nf, nf)
  cor_wsum <- 0
  lambda_share <- numeric(nc)
  components <- vector("list", nc)

  for (ci in seq_len(nc)) {
    block <- blocks[[ci]]
    U <- t4_residual_gene_matrix(block, if (is.null(Zlist)) block$Z else Zlist[[ci]])
    P <- U %*% block$W
    C <- suppressWarnings(stats::cor(P))
    zc_full[, ci] <- t4_atanh(C[cbind(pairs$i, pairs$j)])

    op <- ops[[ci]]
    A <- (P[, op$i, drop = FALSE] - as.matrix(U %*% op$Va)) /
      rep(1 - op$shared_a, each = block$n)
    B <- (P[, op$j, drop = FALSE] - as.matrix(U %*% op$Vb)) /
      rep(1 - op$shared_b, each = block$n)
    zd <- zc_full[, ci]
    zd[op$pair_index] <- t4_atanh(t4_colwise_cor(A, B))
    zd[op$pair_index[!op$eligible]] <- NA_real_
    zc_disj[, ci] <- zd

    # THE GLOBAL RESIDUAL COMPONENTS.
    # Program scores from one bulk library share a large amount of variance that
    # has nothing to do with any specific pair, so most edges are nonzero for one
    # uninteresting reason. The leading components of the program residual space
    # are therefore projected out of BOTH members of every pair and the edges are
    # recomputed, at k = 1 and again at k = 3.
    #
    # TWO PROPERTIES OF THIS ARM THAT CONSTRAIN HOW IT MAY BE READ.
    #  - The components are estimated with the pair included, so a pair that
    #    loads on them loses some of its own signal. A surviving edge is
    #    therefore a conservative claim and never an inflated one.
    #  - Removing k dimensions from p programs forces the leftover covariances to
    #    sum down, so edges that were purely global come back NEGATIVE. Negative
    #    partial edges are an artifact of the projection and are not
    #    interpretable; only positive ones are.
    #  - The components are ESTIMATED, not known. Two programs that are little
    #    more than a readout of the global axis, with no private structure of
    #    their own, both retain the part of it the estimate misses, and they stay
    #    correlated. t4_run_tests.R plants exactly that case: at p = 24 the arm
    #    removes a planted global component from 274 of 275 pairs and the one
    #    survivor is the single pair with no private factor. This is a limit of
    #    the arm, not something a larger k fixes, and it is why the k = 1 count
    #    is never reported without the k = 3 sensitivity beside it.
    Ps <- scale(P)
    Ps[!is.finite(Ps)] <- 0
    sv <- svd(Ps, nu = max(n_components), nv = 0L)
    lambda_share[ci] <- sv$d[[1L]]^2 / sum(sv$d^2)
    components[[ci]] <- sv$u[, seq_len(max(n_components)), drop = FALSE]
    for (ki in seq_along(n_components)) {
      u <- sv$u[, seq_len(n_components[[ki]]), drop = FALSE]
      drop_u <- function(M) M - u %*% crossprod(u, M)
      Cp <- suppressWarnings(stats::cor(drop_u(P)))
      zp <- t4_atanh(Cp[cbind(pairs$i, pairs$j)])
      zp[op$pair_index] <- t4_atanh(t4_colwise_cor(drop_u(A), drop_u(B)))
      zp[op$pair_index[!op$eligible]] <- NA_real_
      zc_part[[ki]][, ci] <- zp
    }

    weight_c[ci] <- block$df - 3L
    Cf <- C; Cf[!is.finite(Cf)] <- 0
    cor_sum <- cor_sum + weight_c[ci] * t4_atanh(Cf)
    cor_wsum <- cor_wsum + weight_c[ci]
  }

  combine <- function(zc) {
    ok <- is.finite(zc)
    w <- matrix(weight_c, np, nc, byrow = TRUE) * ok
    denom <- rowSums(w)
    out <- rowSums(replace(zc, !ok, 0) * w) / denom
    list(z = ifelse(denom > 0, out, NA_real_), weight = denom,
         n_cohorts = rowSums(ok),
         n_same_sign = pmax(rowSums(ok & zc > 0), rowSums(ok & zc < 0)))
  }
  full <- combine(zc_full)
  disj <- combine(zc_disj)
  part <- lapply(zc_part, combine)
  names(part) <- paste0("z_partial_k", n_components)

  Cbar <- tanh(cor_sum / cor_wsum)
  diag(Cbar) <- 1
  eig <- eigen(Cbar, symmetric = TRUE, only.values = TRUE)$values

  agreement <- function(zc) {
    ok <- which(rowSums(is.finite(zc)) == nc)
    if (length(ok) < 3L) return(NA_real_)
    cc <- stats::cor(zc[ok, , drop = FALSE])
    mean(cc[upper.tri(cc)])
  }

  out <- c(
    list(z_full = full$z, z_disjoint = disj$z,
         n_cohorts_disjoint = disj$n_cohorts,
         n_same_sign_disjoint = disj$n_same_sign),
    lapply(part, `[[`, "z"),
    setNames(lapply(part, `[[`, "n_cohorts"), paste0("n_cohorts_k", n_components)),
    setNames(lapply(part, `[[`, "n_same_sign"), paste0("n_same_sign_k", n_components)))
  out$globals <- c(
    mean_abs_z_disjoint = mean(abs(disj$z), na.rm = TRUE),
    mean_abs_z_partial_k1 = mean(abs(part[[1L]]$z), na.rm = TRUE),
    positive_fraction_disjoint = mean(disj$z > 0, na.rm = TRUE),
    positive_fraction_partial_k1 = mean(part[[1L]]$z > 0, na.rm = TRUE),
    lambda1 = eig[[1L]],
    lambda1_share = eig[[1L]] / sum(eig),
    lambda2 = eig[[2L]],
    spectrum_top10_share = sum(eig[seq_len(10L)]) / sum(eig),
    global_component_share = mean(lambda_share),
    cross_cohort_agreement = agreement(zc_disj),
    cross_cohort_agreement_partial_k1 = agreement(zc_part[[1L]]))
  if (want_matrix) {
    out$zc_disjoint <- zc_disj
    out$zc_partial <- zc_part
    out$zc_full <- zc_full
    out$correlation_matrix <- Cbar
    out$eigenvalues <- eig
    out$global_components <- components
  }
  out
}

# Fisher transform with the endpoints pulled off the boundary so a perfectly
# collinear pair does not become infinite.
t4_atanh <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))

# Column-by-column correlation of two matrices with the same shape.
t4_colwise_cor <- function(A, B) {
  Ac <- sweep(A, 2L, colMeans(A), "-")
  Bc <- sweep(B, 2L, colMeans(B), "-")
  num <- colSums(Ac * Bc)
  den <- sqrt(colSums(Ac^2) * colSums(Bc^2))
  ifelse(den > 0, num / den, NA_real_)
}

# ---------------------------------------------------------------------- nulls

# NULL A, the primary. Donors are permuted within cohort independently for each
# gene. Every gene keeps its own residual distribution and every program keeps
# its exact weight vector and its exact shared genes, so the edge correlation
# that survives is precisely the part shared membership alone produces. All
# gene-gene coordination, which is the thing being tested, is destroyed.
t4_permute_genes <- function(Z) {
  g <- nrow(Z); n <- ncol(Z)
  perm <- matrixStats::rowRanks(matrix(stats::runif(g * n), g, n),
                                ties.method = "first")
  out <- matrix(Z[(perm - 1L) * g + seq_len(g)], g, n)
  dimnames(out) <- dimnames(Z)
  out
}

# NULL B, robustness for the leave-shared-out arm. One permutation per gene
# GROUP rather than per gene, so coordination WITHIN a program survives and only
# coordination BETWEEN programs is destroyed. A gene in several programs can
# only take one permutation, so it is assigned to its first program; the null is
# therefore approximate for entangled programs and exact for disjoint ones,
# which is the arm it is used to check.
t4_permute_groups <- function(Z, group) {
  g <- nrow(Z); n <- ncol(Z)
  perms <- matrix(0L, max(group), n)
  for (k in seq_len(max(group))) perms[k, ] <- sample.int(n)
  perm <- perms[group, , drop = FALSE]
  out <- matrix(Z[(perm - 1L) * g + seq_len(g)], g, n)
  dimnames(out) <- dimnames(Z)
  out
}

t4_gene_groups <- function(blocks, weights, features) {
  first <- weights[, .(feature_id = feature_id[[1L]]), by = gene_symbol]
  lapply(blocks, function(block) {
    idx <- match(block$genes, first$gene_symbol)
    grp <- match(first$feature_id[idx], features)
    grp[is.na(grp)] <- 0L
    # genes in no program still need a group; give each its own.
    free <- which(grp == 0L)
    grp[free] <- length(features) + seq_along(free)
    as.integer(grp)
  })
}

# ------------------------------------------------------------ misspecification

# The gate that would have caught the sibling system's retraction. Under a
# residual model saturated in the recorded histology there must be no systematic
# association left between the residual scores and stage.
#
# TWO COLUMNS, BECAUSE THEY MEAN DIFFERENT THINGS.
#  raw_stage_r2 is the R-squared of the UNSCALED residual on the full stage
#  factor. Saturation makes this exactly zero by construction, so a nonzero
#  value is a coding error and nothing else. It is the check that the model
#  really is saturated.
#  stage_r2 uses the leverage-scaled residual that the correlations are actually
#  computed from. Dividing by sqrt(1 - h) is a per-donor rescaling that does not
#  respect the design's column space, so this one is small but not exactly zero.
#  Its reference is the chance level for a 4-df regressor, roughly 4/(n - 1).
t4_stage_leakage <- function(blocks) {
  rows <- lapply(blocks, function(block) {
    Yt <- t(block$Z)
    raw <- Yt - block$X %*% (block$xtx_inv %*% crossprod(block$X, Yt))
    P_raw <- raw %*% block$W
    P <- (raw * block$dinv) %*% block$W
    keep <- which(colSums(is.finite(P)) == nrow(P) & apply(P, 2L, stats::sd) > 0)
    stage <- factor(block$meta$fibrosis_stage)
    r2 <- function(M) suppressWarnings(apply(M, 2L, function(v)
      summary(stats::lm(v ~ stage))$r.squared))
    rho <- suppressWarnings(apply(P[, keep, drop = FALSE], 2L, function(v)
      stats::cor(v, block$meta$fibrosis_stage, method = "spearman")))
    data.table(cohort = block$cohort, feature_id = colnames(P)[keep],
               spearman_stage = rho,
               raw_stage_r2 = r2(P_raw[, keep, drop = FALSE]),
               stage_r2 = r2(P[, keep, drop = FALSE]),
               chance_stage_r2 = (nlevels(stage) - 1L) / (block$n - 1L))
  })
  rbindlist(rows)
}

# ------------------------------------------------------------------- reporting

# Smallest excess correlation, beyond what shared membership produces, that the
# design could have resolved at 80 percent power. Expressed on the correlation
# scale around the pair's own null centre so that a negative on an entangled
# pair is not confused with a negative on a disjoint one.
t4_minimum_detectable_excess <- function(null_mean, null_sd, alpha, power = 0.80) {
  delta <- null_sd * (stats::qnorm(1 - alpha / 2) + stats::qnorm(power))
  tanh(null_mean + delta) - tanh(null_mean)
}

# Louvain partition of the positive excess-coordination graph and its modularity.
# The question it answers is whether the surviving coordination organises into
# blocks or is scattered, and it is compared against the same statistic computed
# on null replicates.
t4_modularity <- function(pairs, weight, features, min_weight = 0) {
  keep <- which(is.finite(weight) & weight > min_weight)
  if (length(keep) < 3L) return(list(modularity = NA_real_, n_communities = NA_integer_,
                                     membership = NULL))
  g <- igraph::graph_from_data_frame(
    data.frame(from = pairs$program_a[keep], to = pairs$program_b[keep],
               weight = weight[keep]),
    directed = FALSE,
    vertices = data.frame(name = features))
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  list(modularity = igraph::modularity(cl), n_communities = length(unique(igraph::membership(cl))),
       n_edges = igraph::ecount(g), graph = g,
       membership = data.table(feature_id = igraph::V(g)$name,
                               community = as.integer(igraph::membership(cl))))
}

# Degree-preserving rewiring null for the modularity above.
#
# The permutation null cannot serve here: under it almost no edge survives BH, so
# the null graph is empty and its modularity is undefined rather than low. The
# question that can be asked is whether the block structure is more than the
# degree sequence implies, which is what rewiring while holding every program's
# degree fixed tests.
t4_rewired_modularity <- function(graph, n_rewirings = 200L) {
  if (is.null(graph)) return(rep(NA_real_, n_rewirings))
  vapply(seq_len(n_rewirings), function(i) {
    g <- igraph::rewire(graph, igraph::keeping_degseq(niter = 10L * igraph::ecount(graph)))
    igraph::modularity(igraph::cluster_louvain(g, weights = igraph::E(g)$weight))
  }, numeric(1))
}
