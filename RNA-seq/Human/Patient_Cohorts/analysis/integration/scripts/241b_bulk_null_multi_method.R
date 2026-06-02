# 241b_bulk_null_multi_method.R
# Phase 1.1b — Re-test bulk F3 sub-state NULL using THREE independent methods
# beyond ConsensusClusterPlus.
#
# Methods:
#   1. PAM (Partitioning Around Medoids) at k=2,3,4 with Spearman distance
#   2. NMF k=2 vs k=1 BIC comparison (NMF::nmf, nrun=20)
#   3. Archetypal analysis with k=3 archetypes (archetypes pkg)
#
# Pre-registered acceptance criteria for "NULL holds" (all three must concur):
#   PAM:        min cluster fraction <15% at every k=2,3,4
#   NMF:        ΔBIC favours k=1 by ≥5 (i.e. BIC_k1 - BIC_k2 ≤ -5; since lower
#               BIC is preferred, we say k=1 is favoured when BIC_k1 < BIC_k2 - 5)
#   Archetypes: posterior mass diffuse (no single archetype >70%, no two >85%)
#
# Inputs:
#   results/integration/merged_dge.rds     (DGEList)
#   results/integration/meta_matched.rds   (sample metadata)
#
# Outputs (to results/granular_staging/):
#   bulk_null_multi_method.csv  — method, metric, value, pass_null, notes
#   bulk_null_multi_method.md   — 1-paragraph verdict

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(matrixStats)
  library(cluster)        # pam, silhouette
  library(NMF)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                          "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT,
                     "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== Phase 1.1b F3 NULL Multi-Method Re-test ==\n")
cat("Output:", OUT_DIR, "\n\n")

# ---------------------------------------------------------------------------
# Ensure archetypes package is available (lazy install if missing).
# Install into the active R lib (writable per CLAUDE).
# ---------------------------------------------------------------------------
have_archetypes <- requireNamespace("archetypes", quietly = TRUE)
if (!have_archetypes) {
  cat("archetypes not installed; attempting CRAN install...\n")
  tryCatch({
    install.packages(
      "archetypes",
      repos = "https://cran.r-project.org",
      quiet = TRUE
    )
    have_archetypes <- requireNamespace("archetypes", quietly = TRUE)
  }, error = function(e) {
    cat("install failed:", conditionMessage(e), "\n")
  })
}
if (!have_archetypes) {
  stop("archetypes package not available; aborting per spec")
}
suppressPackageStartupMessages(library(archetypes))
cat("archetypes loaded.\n\n")

# ---------------------------------------------------------------------------
# Load data and subset to F3
# ---------------------------------------------------------------------------
dge <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta <- meta[match(colnames(dge), meta$sample_id), ]
stopifnot(identical(colnames(dge), meta$sample_id))

f3_idx <- which(!is.na(meta$fibrosis_stage) & meta$fibrosis_stage == 3)
cat(sprintf("F3 samples: %d (across %d cohorts)\n",
            length(f3_idx), length(unique(meta$dataset[f3_idx]))))
print(table(meta$dataset[f3_idx]))

dge_f3 <- dge[, f3_idx]
cohort_vec <- factor(meta$dataset[f3_idx])

# ---------------------------------------------------------------------------
# voom + cohort batch correction + top-2000 variable genes
# (matches Script 241 pooled-analysis pattern exactly)
# ---------------------------------------------------------------------------
keep <- filterByExpr(dge_f3,
                     group = factor(rep("F3", ncol(dge_f3))),
                     min.count = 5)
dge_f3 <- dge_f3[keep, , keep.lib.sizes = FALSE]
dge_f3 <- calcNormFactors(dge_f3)
v <- voom(dge_f3, design = NULL)
expr_raw <- v$E
expr_bc <- removeBatchEffect(expr_raw, batch = cohort_vec)
gene_vars <- rowVars(expr_bc)
top <- order(gene_vars, decreasing = TRUE)[seq_len(min(2000, nrow(expr_bc)))]
expr <- expr_bc[top, ]   # 2000 × 270  (genes × samples)

cat(sprintf("\nExpression matrix (post voom + batch correct + top-var): %d × %d\n",
            nrow(expr), ncol(expr)))

# Sample-by-gene matrix used by some methods
expr_T <- t(expr)        # samples × genes

# Pre-compute Spearman distance once
spearman_dist <- as.dist(1 - cor(expr, method = "spearman"))

# Container for results
results <- list()

# ---------------------------------------------------------------------------
# Method 1: PAM (k=2,3,4) with Spearman distance
# ---------------------------------------------------------------------------
cat("\n--- Method 1: PAM ---\n")
pam_min_frac <- numeric(0)
pam_sil <- numeric(0)
for (k in 2:4) {
  pam_fit <- cluster::pam(spearman_dist, k = k, diss = TRUE)
  cls <- pam_fit$clustering
  tab <- table(cls)
  min_frac <- min(tab) / sum(tab)
  sil <- pam_fit$silinfo$avg.width
  pam_min_frac[as.character(k)] <- min_frac
  pam_sil[as.character(k)] <- sil
  cat(sprintf("k=%d: cluster sizes = %s | min_frac = %.3f | silhouette = %.3f\n",
              k, paste(tab, collapse = "/"), min_frac, sil))
  results[[length(results) + 1]] <- data.frame(
    method = "PAM",
    metric = sprintf("k%d_min_cluster_fraction", k),
    value  = min_frac,
    pass_null = min_frac < 0.15,
    notes  = sprintf("sizes=%s; silhouette=%.3f",
                     paste(tab, collapse = "/"), sil),
    stringsAsFactors = FALSE
  )
  results[[length(results) + 1]] <- data.frame(
    method = "PAM",
    metric = sprintf("k%d_silhouette_avg", k),
    value  = sil,
    pass_null = NA,
    notes  = "diagnostic only",
    stringsAsFactors = FALSE
  )
}
# Aggregate PAM null verdict: NULL holds if ALL k have min_frac < 0.15
pam_pass <- all(pam_min_frac < 0.15)
results[[length(results) + 1]] <- data.frame(
  method = "PAM",
  metric = "AGGREGATE_null_holds",
  value  = max(pam_min_frac),
  pass_null = pam_pass,
  notes  = sprintf("max min_cluster_fraction across k=2,3,4 = %.3f; NULL holds if ALL <0.15",
                   max(pam_min_frac)),
  stringsAsFactors = FALSE
)
cat(sprintf("PAM aggregate: NULL holds = %s\n", pam_pass))

# ---------------------------------------------------------------------------
# Method 2: NMF k=2 vs k=1 BIC
#
# NMF requires non-negative input. voom log-CPM has negatives, so we shift
# by row min and add a small epsilon.
# Strategy:
#   * fit NMF k=2 (nrun=20)
#   * collapse k=2 expression onto a 1-D score (1st PC of the basis-loadings
#     embedding) — equivalent to "rank-1 NMF approximation" in residual space
#   * compute Gaussian BIC for k=1 (single cluster) vs k=2 (two-component
#     Gaussian mixture on the collapsed score)
#   * ΔBIC = BIC_k1 - BIC_k2; ΔBIC ≥ 5 favours k=1
# ---------------------------------------------------------------------------
cat("\n--- Method 2: NMF k=2 vs k=1 BIC ---\n")

# Shift expression to non-negative for NMF input
expr_nn <- expr - min(expr) + 1e-3
cat(sprintf("NMF input matrix: %d × %d (non-negative shift = %.4f)\n",
            nrow(expr_nn), ncol(expr_nn), -min(expr) + 1e-3))

# Fit NMF k=2 with nrun=20 (use brunet default algorithm)
nmf_fit <- tryCatch({
  NMF::nmf(expr_nn, rank = 2, nrun = 20, seed = 42, .options = "v-")
}, error = function(e) {
  cat("NMF::nmf error:", conditionMessage(e), "\n"); NULL
})
stopifnot(!is.null(nmf_fit))

# Collapsed 1-D score: project samples onto difference of two basis loadings.
# coef(nmf_fit) is k × samples; the difference between the two rows captures
# the "axis" between the two NMF components.
H <- coef(nmf_fit)           # 2 × n
score_vec <- as.numeric(H[1, ] - H[2, ])
score_vec <- (score_vec - mean(score_vec)) / sd(score_vec)

# --- Likelihoods on the 1-D score ---
n <- length(score_vec)

# k=1 model: single Gaussian
mu1  <- mean(score_vec)
sd1  <- sd(score_vec)
ll_k1 <- sum(dnorm(score_vec, mean = mu1, sd = sd1, log = TRUE))
p_k1 <- 2                  # mu, sd
bic_k1 <- -2 * ll_k1 + p_k1 * log(n)

# k=2 model: two-component Gaussian mixture (use mixtools::normalmixEM if
# available, else hand-rolled EM).  Light-weight EM here to avoid extra deps.
em2 <- function(x, max_iter = 200, tol = 1e-6) {
  n <- length(x)
  pi  <- c(0.5, 0.5)
  mu  <- quantile(x, c(0.25, 0.75))
  sd2 <- rep(sd(x), 2)
  ll_old <- -Inf
  for (it in seq_len(max_iter)) {
    # E-step: posterior responsibilities
    d1 <- dnorm(x, mu[1], sd2[1])
    d2 <- dnorm(x, mu[2], sd2[2])
    w1 <- pi[1] * d1
    w2 <- pi[2] * d2
    z  <- w1 / (w1 + w2 + 1e-30)
    # M-step
    n1 <- sum(z); n2 <- sum(1 - z)
    pi <- c(n1 / n, n2 / n)
    mu <- c(sum(z * x) / n1, sum((1 - z) * x) / n2)
    sd2 <- c(sqrt(sum(z * (x - mu[1])^2) / n1),
             sqrt(sum((1 - z) * (x - mu[2])^2) / n2))
    sd2 <- pmax(sd2, 1e-6)
    ll <- sum(log(pi[1] * dnorm(x, mu[1], sd2[1]) +
                  pi[2] * dnorm(x, mu[2], sd2[2]) + 1e-30))
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(pi = pi, mu = mu, sd = sd2, loglik = ll, iter = it)
}

em_fit <- em2(score_vec)
ll_k2 <- em_fit$loglik
p_k2 <- 5                  # 2 mu + 2 sd + 1 mixing weight (other constrained)
bic_k2 <- -2 * ll_k2 + p_k2 * log(n)
delta_bic <- bic_k1 - bic_k2    # >0 means k=2 has lower BIC

cat(sprintf("BIC k=1: %.3f | BIC k=2: %.3f | ΔBIC (k1-k2): %.3f\n",
            bic_k1, bic_k2, delta_bic))
nmf_pass <- delta_bic <= -5     # k=1 favoured when delta_bic <= -5
                                # (i.e. BIC_k1 lower than BIC_k2 by ≥5)
cat(sprintf("NMF NULL holds (k=1 favoured by ≥5 BIC): %s\n", nmf_pass))

results[[length(results) + 1]] <- data.frame(
  method = "NMF",
  metric = "BIC_k1",
  value  = bic_k1,
  pass_null = NA,
  notes  = sprintf("Gaussian on 1D collapsed score n=%d", n),
  stringsAsFactors = FALSE
)
results[[length(results) + 1]] <- data.frame(
  method = "NMF",
  metric = "BIC_k2",
  value  = bic_k2,
  pass_null = NA,
  notes  = "2-component GMM via EM",
  stringsAsFactors = FALSE
)
results[[length(results) + 1]] <- data.frame(
  method = "NMF",
  metric = "delta_BIC_k1_minus_k2",
  value  = delta_bic,
  pass_null = nmf_pass,
  notes  = sprintf("k=1 favoured when ≤-5; observed %.3f", delta_bic),
  stringsAsFactors = FALSE
)
results[[length(results) + 1]] <- data.frame(
  method = "NMF",
  metric = "AGGREGATE_null_holds",
  value  = delta_bic,
  pass_null = nmf_pass,
  notes  = "ΔBIC ≤ -5 ⇒ NULL holds",
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Method 3: Archetypal analysis (k=3)
#
# Archetypes work in sample-space; pass samples × top-genes (n × p). Each
# sample gets a posterior mixture (alpha) over 3 archetypes. NULL holds iff
# mass is diffuse: no single archetype dominates >70% of mass and no two
# archetypes >85%.
# Diagnostic: per-sample max alpha and mean per-archetype mass share.
# ---------------------------------------------------------------------------
cat("\n--- Method 3: Archetypes (k=3) ---\n")
arch_in <- expr_T            # samples × genes (270 × 2000)
arch_fit <- tryCatch({
  archetypes::archetypes(arch_in, k = 3, verbose = FALSE)
}, error = function(e) {
  cat("archetypes error:", conditionMessage(e), "\n"); NULL
})
stopifnot(!is.null(arch_fit))

alpha <- arch_fit$alphas      # n × 3
stopifnot(ncol(alpha) == 3)
arch_mass <- colMeans(alpha)  # mean fraction of mass each archetype carries
arch_mass_sorted <- sort(arch_mass, decreasing = TRUE)
top1 <- arch_mass_sorted[1]
top2 <- arch_mass_sorted[1] + arch_mass_sorted[2]
arch_rss <- arch_fit$rss

cat(sprintf("Archetype mean alpha mass (sorted): %s\n",
            paste(sprintf("%.3f", arch_mass_sorted), collapse = " | ")))
cat(sprintf("Top-1 mass: %.3f | Top-2 mass: %.3f | RSS: %.4f\n",
            top1, top2, arch_rss))
arch_pass <- (top1 < 0.70) && (top2 < 0.85)
cat(sprintf("Archetypes NULL holds (diffuse mass): %s\n", arch_pass))

results[[length(results) + 1]] <- data.frame(
  method = "Archetypes",
  metric = "top1_mean_alpha",
  value  = top1,
  pass_null = (top1 < 0.70),
  notes  = "max archetype mean mass; <0.70 ⇒ no single dominant archetype",
  stringsAsFactors = FALSE
)
results[[length(results) + 1]] <- data.frame(
  method = "Archetypes",
  metric = "top2_mean_alpha",
  value  = top2,
  pass_null = (top2 < 0.85),
  notes  = "sum of top-2 archetype mean masses; <0.85 ⇒ third matters",
  stringsAsFactors = FALSE
)
results[[length(results) + 1]] <- data.frame(
  method = "Archetypes",
  metric = "rss",
  value  = arch_rss,
  pass_null = NA,
  notes  = "diagnostic: residual sum of squares of fit",
  stringsAsFactors = FALSE
)
results[[length(results) + 1]] <- data.frame(
  method = "Archetypes",
  metric = "archetype_mass_3rd",
  value  = arch_mass_sorted[3],
  pass_null = NA,
  notes  = "smallest mean alpha (>0.15 supports diffuse mass)",
  stringsAsFactors = FALSE
)
results[[length(results) + 1]] <- data.frame(
  method = "Archetypes",
  metric = "AGGREGATE_null_holds",
  value  = top1,
  pass_null = arch_pass,
  notes  = sprintf("top1=%.3f, top2=%.3f; NULL = (top1<0.70 AND top2<0.85)",
                   top1, top2),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Combine + verdict
# ---------------------------------------------------------------------------
out_df <- do.call(rbind, results)
csv_path <- file.path(OUT_DIR, "bulk_null_multi_method.csv")
write.csv(out_df, csv_path, row.names = FALSE)
cat(sprintf("\nWrote %s (%d rows)\n", csv_path, nrow(out_df)))

# Verdict logic
all_pass  <- pam_pass && nmf_pass && arch_pass
all_fail  <- (!pam_pass) && (!nmf_pass) && (!arch_pass)
verdict <- if (all_pass) {
  "NULL HOLDS"
} else if (all_fail) {
  "NULL COLLAPSES"
} else {
  "METHOD-DEPENDENT"
}

# 1-paragraph md summary
md_path <- file.path(OUT_DIR, "bulk_null_multi_method.md")
md_lines <- c(
  "# F3 sub-state bulk NULL — three-method re-test (Phase 1.1b)",
  "",
  sprintf("**Verdict: %s**", verdict),
  "",
  sprintf(paste0(
    "Re-tested the F3 sub-state NULL on n=%d F3 samples ",
    "(top-2000 variable genes, voom + cohort-batch corrected) using three ",
    "independent methods beyond ConsensusClusterPlus. (1) PAM with ",
    "Spearman distance: min cluster fraction was %.3f / %.3f / %.3f at ",
    "k=2/3/4 (NULL gate <0.15 ⇒ %s). (2) NMF k=2 fitted with nrun=20; ",
    "the 1-D collapsed score gave BIC_k1=%.2f vs BIC_k2=%.2f, ΔBIC=%.2f ",
    "(NULL gate ≤-5 ⇒ %s). (3) Archetypal analysis k=3: mean posterior ",
    "mass per archetype was %s (sorted), with top-1=%.3f and top-2=%.3f ",
    "(NULL gate top-1<0.70 AND top-2<0.85 ⇒ %s). Verdict: **%s** ",
    "(PAM=%s, NMF=%s, Archetypes=%s)."
  ),
    ncol(expr),
    pam_min_frac["2"], pam_min_frac["3"], pam_min_frac["4"],
    pam_pass,
    bic_k1, bic_k2, delta_bic, nmf_pass,
    paste(sprintf("%.3f", arch_mass_sorted), collapse = "/"),
    top1, top2, arch_pass,
    verdict, pam_pass, nmf_pass, arch_pass)
)
writeLines(md_lines, md_path)
cat(sprintf("Wrote %s\n", md_path))

cat(sprintf("\n=== VERDICT: %s ===\n", verdict))
cat(sprintf("PAM pass=%s | NMF pass=%s | Archetypes pass=%s\n",
            pam_pass, nmf_pass, arch_pass))
cat("\nDone.\n")
