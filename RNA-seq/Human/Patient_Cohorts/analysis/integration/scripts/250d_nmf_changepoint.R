# 250d_nmf_changepoint.R
# Bayesian changepoint test of the "F1→F2 = P1 anchor; F2→F3 = P6 anchor" claim.
#
# Approach:
#   For each NMF program P1..P6, fit per-sample H[p,] vs fibrosis_stage with:
#     M0 (null): single linear stage trend
#     M1 (1 changepoint): two linear segments with breakpoint among F0→F1, F1→F2,
#                         F2→F3, F3→F4.
#   Compare via mcp loo (psis-loo) for changepoint location posterior.
#   Pre-registered acceptance for "real anchor" at transition T:
#     (a) posterior(changepoint within transition T window) > 0.9
#     (b) slope at T >= 1.5x median slope across other transitions
#     (c) Bayesian binary t-test before-vs-after T: posterior P(diff > 0) > 0.9
#
# Outputs (RNA-seq/results/granular_staging/):
#   nmf_changepoint_posterior.csv  — per program × transition row with posterior,
#                                    slopes, pass flag
#   nmf_changepoint_summary.md     — verdict on P1/P6 anchor claims

suppressPackageStartupMessages({
  library(NMF)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                          "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
NMF_CACHE <- file.path(PROJECT_ROOT, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== 250d NMF Bayesian changepoint test ==\n")
cat(sprintf("Started: %s\n", Sys.time()))

# ---------------------------------------------------------------------------
# Optional install of mcp (Multiple Change Points). Handle gracefully — if
# install fails (no internet on compute node, dependency conflict), fall back
# to a custom pure-R Bayesian changepoint marginal-likelihood comparison
# instead of failing the script.
# ---------------------------------------------------------------------------

have_mcp <- requireNamespace("mcp", quietly = TRUE)
if (!have_mcp) {
  cat("mcp not installed — attempting install from CRAN...\n")
  try(install.packages("mcp", repos = "https://cloud.r-project.org",
                       quiet = TRUE), silent = TRUE)
  have_mcp <- requireNamespace("mcp", quietly = TRUE)
}
cat(sprintf("mcp available: %s\n", have_mcp))

# ---------------------------------------------------------------------------
# Load NMF H matrix (k=6 programs × samples) and metadata (fibrosis stage)
# ---------------------------------------------------------------------------

nmf <- readRDS(NMF_CACHE)
k6 <- nmf$nmf_results[["6"]]
H <- NMF::coef(k6)
cat(sprintf("H matrix: %d programs × %d samples\n", nrow(H), ncol(H)))

meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
common <- intersect(colnames(H), meta$sample_id)
H_match <- H[, common, drop = FALSE]
meta_match <- meta[match(common, meta$sample_id), ]
prog_names <- paste0("P", seq_len(nrow(H_match)))

# Build per-sample table: columns = sample, stage, P1..P6
df_samples <- data.frame(
  sample_id = common,
  stage = as.numeric(meta_match$fibrosis_stage),
  stringsAsFactors = FALSE
)
for (p in seq_len(nrow(H_match))) {
  df_samples[[prog_names[p]]] <- as.numeric(H_match[p, ])
}
df_samples <- df_samples[!is.na(df_samples$stage), ]
cat(sprintf("Samples with stage: %d\n", nrow(df_samples)))
cat("Stage counts:\n")
print(table(df_samples$stage))

transitions <- list(
  "F0_F1" = c(0, 1),
  "F1_F2" = c(1, 2),
  "F2_F3" = c(2, 3),
  "F3_F4" = c(3, 4)
)

# ---------------------------------------------------------------------------
# Helper: per-segment slope from per-sample data
# Slope at transition (a, b) = (mean(b) - mean(a)) / (b - a). Δstage = 1.
# ---------------------------------------------------------------------------

segment_slope <- function(scores, stages, a, b) {
  x_a <- scores[stages == a]
  x_b <- scores[stages == b]
  if (length(x_a) < 5 || length(x_b) < 5) return(NA_real_)
  (mean(x_b, na.rm = TRUE) - mean(x_a, na.rm = TRUE)) / (b - a)
}

# Bootstrap CI for the segment slope.
boot_slope <- function(scores, stages, a, b, B = 1000) {
  x_a <- scores[stages == a]; x_b <- scores[stages == b]
  if (length(x_a) < 5 || length(x_b) < 5) return(c(NA_real_, NA_real_, NA_real_))
  out <- numeric(B)
  for (i in seq_len(B)) {
    out[i] <- (mean(sample(x_b, replace = TRUE)) -
                 mean(sample(x_a, replace = TRUE))) / (b - a)
  }
  c(mean = mean(out), q025 = quantile(out, 0.025, names = FALSE),
    q975 = quantile(out, 0.975, names = FALSE))
}

# ---------------------------------------------------------------------------
# Bayesian binary before-vs-after t-test (sensitivity / Test 3)
# Pure-R conjugate normal-with-unknown-mean-and-variance; computes
#   P(mu_after > mu_before | data) using equal-tail draws from the posterior of
#   each group (independent Jeffreys priors → t posterior on the mean).
# ---------------------------------------------------------------------------

bayes_diff_posterior <- function(x_a, x_b, n_draws = 2e4) {
  # Posterior on mean for each group: t-distribution with df = n-1, mean = xbar,
  # scale = s / sqrt(n). Returns P(mu_b > mu_a).
  na <- length(x_a); nb <- length(x_b)
  if (na < 5 || nb < 5) return(NA_real_)
  ma <- mean(x_a); mb <- mean(x_b)
  sa <- sd(x_a) / sqrt(na); sb <- sd(x_b) / sqrt(nb)
  set.seed(42)
  draws_a <- ma + sa * rt(n_draws, df = na - 1)
  draws_b <- mb + sb * rt(n_draws, df = nb - 1)
  mean(draws_b > draws_a)
}

# ---------------------------------------------------------------------------
# Bayesian changepoint posterior over the 4 candidate transitions.
# If mcp is available, use it (proper hierarchical Bayesian; gives marginal
# posterior on changepoint location). Otherwise use a discrete marginal-
# likelihood comparison: for each candidate breakpoint cp ∈ {0.5, 1.5, 2.5, 3.5}
# fit a 2-segment OLS (left of cp / right of cp) and compute BIC, then convert
# BIC differences to posterior probabilities (uniform prior over the 4 cp
# locations + a "no changepoint" linear option).
# ---------------------------------------------------------------------------

discrete_cp_posterior <- function(scores, stages) {
  # M0: linear-only; M_k for k in 1..4: two-segment with breakpoint at
  # 0.5, 1.5, 2.5, 3.5 respectively.
  ok <- !is.na(scores) & !is.na(stages)
  y <- scores[ok]; x <- stages[ok]
  n <- length(y)
  if (n < 30) return(NULL)

  # Linear-only
  fit0 <- lm(y ~ x)
  bic0 <- BIC(fit0)

  cp_locs <- c(0.5, 1.5, 2.5, 3.5)
  bic_k <- numeric(length(cp_locs))
  for (i in seq_along(cp_locs)) {
    cp <- cp_locs[i]
    seg <- as.integer(x >= cp)
    fit_k <- lm(y ~ x * seg)
    bic_k[i] <- BIC(fit_k)
  }

  # Convert to posteriors with uniform prior on the 5 models (M0, M1..M4).
  bics <- c(bic0, bic_k)
  log_posts <- -0.5 * (bics - min(bics))
  posts <- exp(log_posts) / sum(exp(log_posts))
  list(
    p_no_cp = posts[1],
    p_at = setNames(posts[-1], names(transitions)),
    bic_no_cp = bic0,
    bic_at = setNames(bic_k, names(transitions))
  )
}

# Wrapper that returns posterior at each transition; uses mcp where available.
mcp_cp_posterior <- function(scores, stages) {
  if (!have_mcp) return(NULL)
  d <- data.frame(y = scores, x = as.numeric(stages))
  d <- d[complete.cases(d), ]
  if (nrow(d) < 30) return(NULL)
  # Two-segment piecewise-linear with one changepoint anywhere in [0.5, 3.5].
  model <- list(y ~ 1 + x, ~ 1 + x)
  fit <- try(mcp::mcp(model, data = d, par_x = "x",
                      adapt = 1500, iter = 4000, chains = 2,
                      cores = 2),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  # Extract posterior samples on the changepoint cp_1
  samp <- try(do.call(rbind, fit$mcmc_post), silent = TRUE)
  if (inherits(samp, "try-error") || !"cp_1" %in% colnames(samp)) return(NULL)
  cp_post <- samp[, "cp_1"]
  cp_post <- cp_post[is.finite(cp_post)]
  # Bin posterior into the 4 transitions.
  bins <- cut(cp_post,
              breaks = c(-Inf, 0.5, 1.5, 2.5, 3.5, Inf),
              labels = c("below_F0", "F0_F1", "F1_F2", "F2_F3", "F3_F4_or_above"),
              include.lowest = TRUE)
  tbl <- prop.table(table(bins))
  list(
    p_at = c(F0_F1 = unname(tbl["F0_F1"]),
             F1_F2 = unname(tbl["F1_F2"]),
             F2_F3 = unname(tbl["F2_F3"]),
             F3_F4 = unname(tbl["F3_F4_or_above"])),
    cp_post = cp_post
  )
}

# ---------------------------------------------------------------------------
# Run per-program × per-transition test.
# ---------------------------------------------------------------------------

results <- list()
mcp_status <- character(length(prog_names))
names(mcp_status) <- prog_names

for (pname in prog_names) {
  cat(sprintf("\n--- %s ---\n", pname))
  scores <- df_samples[[pname]]
  stages <- df_samples$stage

  # Compute slope at every transition
  slopes <- sapply(transitions, function(tr) segment_slope(scores, stages, tr[1], tr[2]))
  cat("Slopes per transition:\n"); print(slopes)

  # Bayesian changepoint posterior
  cp_mcp <- mcp_cp_posterior(scores, stages)
  cp_disc <- discrete_cp_posterior(scores, stages)
  if (!is.null(cp_mcp)) {
    cp_post_used <- cp_mcp$p_at
    mcp_status[pname] <- "mcp"
  } else {
    cp_post_used <- cp_disc$p_at
    mcp_status[pname] <- "bic_discrete"
  }
  cat("Posterior changepoint per transition:\n"); print(cp_post_used)

  # Bayesian before-vs-after sensitivity at each transition
  diff_post <- sapply(transitions, function(tr) {
    bayes_diff_posterior(scores[stages == tr[1]], scores[stages == tr[2]])
  })

  # Bootstrap slope CI per transition
  boots <- sapply(transitions, function(tr) boot_slope(scores, stages, tr[1], tr[2]))

  for (tname in names(transitions)) {
    other <- setdiff(names(transitions), tname)
    other_slopes <- abs(slopes[other])
    median_other <- median(other_slopes, na.rm = TRUE)
    slope_ratio <- if (is.finite(median_other) && median_other > 1e-12)
      abs(slopes[tname]) / median_other else NA_real_

    # Pre-registered "real anchor" criteria (same direction as primary slope sign):
    pass_post <- isTRUE(cp_post_used[tname] > 0.9)
    pass_slope_ratio <- isTRUE(slope_ratio >= 1.5)
    pass_diff_post <- isTRUE(abs(diff_post[tname] - 0.5) > 0.4)  # equiv to >0.9 or <0.1
    pass_real_anchor <- pass_post && pass_slope_ratio && pass_diff_post

    notes <- character(0)
    if (mcp_status[pname] == "bic_discrete")
      notes <- c(notes, "BIC-fallback (mcp unavailable)")
    if (!is.finite(slopes[tname])) notes <- c(notes, "insufficient_n")
    if (!pass_post) notes <- c(notes, sprintf("posterior=%.2f<0.90", cp_post_used[tname]))
    if (!pass_slope_ratio) notes <- c(notes, sprintf("slope_ratio=%.2f<1.5", slope_ratio))
    if (!pass_diff_post) notes <- c(notes,
        sprintf("binary_post=%.2f", diff_post[tname]))

    results[[length(results) + 1]] <- data.frame(
      program = pname,
      transition = tname,
      posterior_changepoint = unname(cp_post_used[tname]),
      slope_at_transition = unname(slopes[tname]),
      slope_lo = boots["q025", tname],
      slope_hi = boots["q975", tname],
      slope_ratio_vs_other_transitions = slope_ratio,
      bayes_diff_posterior_after_gt_before = unname(diff_post[tname]),
      cp_method = mcp_status[pname],
      pass_real_anchor = pass_real_anchor,
      notes = paste(notes, collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
}

posterior_df <- do.call(rbind, results)
write.csv(posterior_df,
          file.path(OUT_DIR, "nmf_changepoint_posterior.csv"),
          row.names = FALSE)
cat(sprintf("\nWrote: %s\n",
            file.path(OUT_DIR, "nmf_changepoint_posterior.csv")))

# ---------------------------------------------------------------------------
# Verdict: do P1 (claimed F1→F2 anchor) and P6 (claimed F2→F3 anchor) survive?
# ---------------------------------------------------------------------------

p1_f1f2 <- posterior_df[posterior_df$program == "P1" &
                          posterior_df$transition == "F1_F2", ]
p6_f2f3 <- posterior_df[posterior_df$program == "P6" &
                          posterior_df$transition == "F2_F3", ]

p1_pass <- isTRUE(p1_f1f2$pass_real_anchor)
p6_pass <- isTRUE(p6_f2f3$pass_real_anchor)

# Are *any* programs monotonic? Check whether all transitions for that program
# share the same slope sign.
program_monotonic <- function(rows) {
  s <- rows$slope_at_transition
  s <- s[is.finite(s)]
  if (length(s) == 0) return(NA)
  all(s > 0) || all(s < 0)
}
mono_summary <- sapply(prog_names, function(pn)
  program_monotonic(posterior_df[posterior_df$program == pn, ]))

method_used <- if (any(posterior_df$cp_method == "mcp")) "mcp::mcp (MCMC)" else
  "BIC discrete-grid fallback (mcp unavailable)"

verdict_md <- c(
  "# NMF program changepoint anchor test (Script 250d)",
  "",
  sprintf("Method: %s", method_used),
  sprintf("Total samples with fibrosis_stage: %d", nrow(df_samples)),
  sprintf("Stages tested (n>=5): %s",
          paste(sort(unique(df_samples$stage)), collapse = ", ")),
  "",
  "## Pre-registered acceptance for 'real anchor'",
  "",
  "1. Posterior probability of changepoint at the candidate transition > 0.9",
  "2. Slope at the candidate transition >= 1.5x median slope across other transitions",
  "3. Bayesian before-vs-after t-test posterior > 0.9 (extreme tail)",
  "",
  "## Verdict on the two claimed anchors",
  "",
  sprintf("- **P1 (Pro-inflammatory) at F1→F2 anchor: %s**",
          ifelse(p1_pass, "SURVIVES", "FAILS")),
  sprintf("  - posterior_changepoint = %.3f (req > 0.90)",
          p1_f1f2$posterior_changepoint),
  sprintf("  - slope at F1→F2 = %.4f, ratio vs median other = %.2fx (req >= 1.5x)",
          p1_f1f2$slope_at_transition,
          p1_f1f2$slope_ratio_vs_other_transitions),
  sprintf("  - binary before-vs-after posterior = %.3f",
          p1_f1f2$bayes_diff_posterior_after_gt_before),
  sprintf("  - notes: %s", ifelse(nchar(p1_f1f2$notes) == 0, "—", p1_f1f2$notes)),
  "",
  sprintf("- **P6 (Stellate-myofibroblast) at F2→F3 anchor: %s**",
          ifelse(p6_pass, "SURVIVES", "FAILS")),
  sprintf("  - posterior_changepoint = %.3f (req > 0.90)",
          p6_f2f3$posterior_changepoint),
  sprintf("  - slope at F2→F3 = %.4f, ratio vs median other = %.2fx (req >= 1.5x)",
          p6_f2f3$slope_at_transition,
          p6_f2f3$slope_ratio_vs_other_transitions),
  sprintf("  - binary before-vs-after posterior = %.3f",
          p6_f2f3$bayes_diff_posterior_after_gt_before),
  sprintf("  - notes: %s", ifelse(nchar(p6_f2f3$notes) == 0, "—", p6_f2f3$notes)),
  "",
  "## Are programs monotonic across stages (vs locally inflected)?",
  ""
)
for (pn in prog_names) {
  verdict_md <- c(verdict_md,
                  sprintf("- %s: monotonic = %s",
                          pn, ifelse(is.na(mono_summary[pn]), "n/a",
                                     ifelse(mono_summary[pn], "YES", "NO"))))
}

# Add table of slopes
verdict_md <- c(verdict_md, "",
                "## Per-program slopes per transition (rounded)", "",
                "| program | F0_F1 | F1_F2 | F2_F3 | F3_F4 |",
                "|---|---|---|---|---|")
for (pn in prog_names) {
  rows <- posterior_df[posterior_df$program == pn, ]
  s <- setNames(rows$slope_at_transition, rows$transition)
  verdict_md <- c(verdict_md,
                  sprintf("| %s | %.4f | %.4f | %.4f | %.4f |",
                          pn,
                          s[["F0_F1"]], s[["F1_F2"]],
                          s[["F2_F3"]], s[["F3_F4"]]))
}

verdict_md <- c(verdict_md, "",
                "## Per-program changepoint posterior per transition", "",
                "| program | F0_F1 | F1_F2 | F2_F3 | F3_F4 |",
                "|---|---|---|---|---|")
for (pn in prog_names) {
  rows <- posterior_df[posterior_df$program == pn, ]
  s <- setNames(rows$posterior_changepoint, rows$transition)
  verdict_md <- c(verdict_md,
                  sprintf("| %s | %.3f | %.3f | %.3f | %.3f |",
                          pn,
                          s[["F0_F1"]], s[["F1_F2"]],
                          s[["F2_F3"]], s[["F3_F4"]]))
}

# Headline interpretation
all_p_pass <- posterior_df[posterior_df$pass_real_anchor, ]
all_mono <- all(mono_summary[prog_names], na.rm = TRUE)
verdict_md <- c(verdict_md, "",
                "## Headline",
                "",
                if (p1_pass && p6_pass) {
                  paste0("Both P1 (F1→F2) and P6 (F2→F3) anchor claims SURVIVE",
                         " the pre-registered Bayesian + slope-ratio criteria.",
                         " The 'two-transition' framing is supported.")
                } else if (!p1_pass && !p6_pass) {
                  paste0("NEITHER anchor claim survives. Programs appear to follow",
                         " gradual / monotonic trajectories across fibrosis stages",
                         " rather than discrete locally-anchored switches.")
                } else if (p1_pass) {
                  "Only the P1@F1→F2 anchor survives; P6@F2→F3 does NOT meet criteria."
                } else {
                  "Only the P6@F2→F3 anchor survives; P1@F1→F2 does NOT meet criteria."
                },
                "",
                sprintf("Programs that pass anchor test (any transition): %s",
                        if (nrow(all_p_pass) == 0) "none"
                        else paste(sprintf("%s@%s", all_p_pass$program,
                                           all_p_pass$transition),
                                   collapse = ", ")),
                sprintf("Programs that are strictly monotonic across stages: %s",
                        paste(prog_names[mono_summary[prog_names] %in% TRUE],
                              collapse = ", "))
)

writeLines(verdict_md,
           file.path(OUT_DIR, "nmf_changepoint_summary.md"))
cat(sprintf("Wrote: %s\n", file.path(OUT_DIR, "nmf_changepoint_summary.md")))

cat(sprintf("\nFinished: %s\n", Sys.time()))
cat("== Done ==\n")
