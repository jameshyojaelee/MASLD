#!/usr/bin/env Rscript
# sex_v3/11b_perm_aggregate.R
# ---------------------------------------------------------------------------
# Pillar 2 aggregator — combine per-rep permutation CSVs into a permutation
# FDR table. Computes:
#   • p_emp_per_gene = (1 + Σ I(|t_perm| ≥ |t_obs|)) / (1 + n_reps)
#   • p_emp_pooled   = (1 + Σ I(|t_perm| ≥ |t_obs|)) / (1 + n_reps × n_genes)
#   • q_emp_per_gene = Storey q on p_emp_per_gene
#   • q_emp_pooled   = Storey q on p_emp_pooled
#
# Inputs: intermediates/perm_v6/rep_*.csv + interaction_classifier_v5_random.csv
# Output: perm_fdr_v6.csv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(qvalue) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")
IDIR  <- file.path(SEXV3, "intermediates")

# Contrast routing (sex_v3_utils.R::contrast_paths) — overrides SEXV3 / IDIR
# when CONTRAST_NAME != "disease_vs_ctrl"; default preserves legacy layout.
if (!exists("contrast_paths")) source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT",
             "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3/sex_v3_utils.R"))
.cpaths <- contrast_paths()
SEXV3 <- .cpaths$sexv3
IDIR  <- .cpaths$idir
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)
PERM_DIR <- file.path(IDIR, "perm_v6")
OUT_CSV  <- file.path(SEXV3, "perm_fdr_v6.csv")
RANDOM_CSV <- file.path(SEXV3, "interaction_classifier_v5_random.csv")
V5_CSV     <- file.path(SEXV3, "interaction_classifier_v5.csv")

cat("=== sex_v6 P2 aggregator ===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# ---- Load observed |t_int| from P1 random-slope refit ----
if (!file.exists(RANDOM_CSV)) stop("P1 output missing: ", RANDOM_CSV)
obs <- fread(RANDOM_CSV, select = c("gene", "t_int_random", "p_int_random",
                                     "padj_int_5k_random"))
setnames(obs, c("gene", "t_int_obs", "p_dream_obs", "padj_int_5k_random"))
obs[, abs_t_obs := abs(t_int_obs)]
cat("Observed genes:", nrow(obs), "\n")

# ---- Read per-rep CSVs ----
rep_files <- list.files(PERM_DIR, pattern = "^rep_\\d+\\.csv$", full.names = TRUE)
n_files <- length(rep_files)
cat("rep_*.csv found:", n_files, "\n")
stopifnot(n_files >= 1)

# Pooled-null tail count requires the full t_perm vector across all reps × genes.
# We use streaming counts to avoid blowing memory.
gene_count <- setNames(integer(nrow(obs)), obs$gene)
abs_t_lookup <- setNames(obs$abs_t_obs, obs$gene)
n_reps_used <- 0L
all_t_perm <- vector("list", n_files)  # for pooled-null

for (i in seq_along(rep_files)) {
  f <- rep_files[i]
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || all(is.na(d$gene))) next
  if (!all(c("gene", "t_int_perm") %in% names(d))) next
  d <- d[!is.na(t_int_perm) & !is.na(gene)]
  if (nrow(d) == 0) next
  # Per-gene exceedance count
  abs_t_perm <- abs(d$t_int_perm)
  abs_t_obs_aligned <- abs_t_lookup[d$gene]
  ok <- !is.na(abs_t_obs_aligned)
  exceed <- abs_t_perm[ok] >= abs_t_obs_aligned[ok]
  hit_genes <- d$gene[ok][exceed]
  if (length(hit_genes) > 0) {
    tbl <- table(hit_genes)
    gene_count[names(tbl)] <- gene_count[names(tbl)] + as.integer(tbl)
  }
  all_t_perm[[i]] <- abs_t_perm
  n_reps_used <- n_reps_used + 1L
}
cat("Valid reps:", n_reps_used, "/", n_files, "\n")
stopifnot(n_reps_used >= 1)

# ---- p_emp_per_gene ----
p_emp_per_gene <- (1 + gene_count) / (1 + n_reps_used)
p_emp_per_gene <- pmin(p_emp_per_gene, 1)

# ---- p_emp_pooled (per gene: P(|t_perm| ≥ |t_obs|) under pooled null) ----
pool <- unlist(all_t_perm, use.names = FALSE)
pool <- pool[!is.na(pool)]
pool_sorted <- sort(pool, decreasing = TRUE)
n_pool <- length(pool_sorted)
cat("Pooled null t values:", n_pool, "\n")
# For each abs_t_obs, count how many pool values are >= abs_t_obs
# Use findInterval on sorted-decreasing → equivalent to sum(pool >= x)
# pool_sorted descending; rank from a sorted-ascending search is easier:
pool_asc <- sort(pool)  # ascending
abs_t_obs <- obs$abs_t_obs
# count of pool >= x = n_pool - (rank of x in pool_asc with right-search)
idx <- findInterval(abs_t_obs, pool_asc, left.open = TRUE)
n_ge <- n_pool - idx
p_emp_pooled <- (1 + n_ge) / (1 + n_pool)
p_emp_pooled <- pmin(p_emp_pooled, 1)

# ---- Storey q-values ----
qv_safe <- function(p) {
  ok <- !is.na(p) & p >= 0 & p <= 1
  out <- rep(NA_real_, length(p))
  if (sum(ok) < 10) return(out)
  q <- tryCatch(qvalue::qvalue(p[ok], pi0 = NULL)$qvalues,
                error = function(e) tryCatch(
                  qvalue::qvalue(p[ok], pi0 = 1)$qvalues,
                  error = function(e2) p.adjust(p[ok], method = "BH")))
  out[ok] <- q
  out
}
q_emp_per_gene <- qv_safe(p_emp_per_gene)
q_emp_pooled   <- qv_safe(p_emp_pooled)

obs[, p_emp_per_gene := p_emp_per_gene]
obs[, p_emp_pooled   := p_emp_pooled]
obs[, q_emp_per_gene := q_emp_per_gene]
obs[, q_emp_pooled   := q_emp_pooled]
obs[, n_perm_used    := n_reps_used]
obs[, n_exceed       := gene_count[obs$gene]]

# ---- agree_dream_vs_perm: do parametric p_dream and empirical p_emp rank-correlate? ----
ok_corr <- !is.na(obs$p_dream_obs) & !is.na(obs$p_emp_per_gene)
spearman_full <- cor(obs$p_dream_obs[ok_corr], obs$p_emp_per_gene[ok_corr],
                     method = "spearman")
# Tier-1 universe: those with padj_int_5k_random not NA (in-Tier-1)
in_tier1 <- !is.na(obs$padj_int_5k_random)
spearman_t1 <- cor(obs$p_dream_obs[in_tier1 & ok_corr],
                   obs$p_emp_per_gene[in_tier1 & ok_corr],
                   method = "spearman")
cat("\nSpearman(p_dream, p_emp_per_gene):\n")
cat("  full:", round(spearman_full, 3), "\n")
cat("  Tier-1:", round(spearman_t1, 3),
    " (n =", sum(in_tier1 & ok_corr), ")\n")
acceptance_pass <- !is.na(spearman_t1) && spearman_t1 >= 0.90
cat("Acceptance gate (Tier-1 Spearman ≥ 0.90):",
    if (acceptance_pass) "PASS" else "FAIL", "\n")
obs[, agree_dream_vs_perm := spearman_t1]

# ---- class_v5_perm — re-run v5 decision tree using q_emp_per_gene where padj_int_5k was ----
# We attach the v5 decision tree (using per-gene q_emp_per_gene as Tier-1 "padj_int_5k") via
# the existing interaction_classifier_v5_random.csv columns.
random <- fread(RANDOM_CSV, select = c("gene", "beta_F", "beta_M", "beta_F_strat",
                                        "beta_M_strat", "lfsr_F_strat", "lfsr_M_strat",
                                        "power_M_at_F", "chr"))
obs2 <- merge(obs, random, by = "gene", all.x = TRUE)
classify_v5 <- function(bF, bM, pI, lF, lM, pow, ch) {
  if (is.na(pI) || is.na(bF) || is.na(bM)) {
    if (!is.na(lF) && !is.na(lM)) {
      if (lF < 0.05 && lM < 0.05 && sign(bF) == sign(bM)) return("Concordant")
      if (lF < 0.05 && lM < 0.05 && sign(bF) != sign(bM)) return("Divergent")
    }
    return("Not_DEG_or_NotInUniverse")
  }
  same_sign <- (sign(bF) == sign(bM)) && bF != 0 && bM != 0
  opp_sign  <- (sign(bF) != sign(bM)) && bF != 0 && bM != 0
  if (pI >= 0.20) {
    if (!is.na(lF) && !is.na(lM)) {
      if (lF < 0.05 && lM < 0.05 && same_sign) return("Concordant")
      if (lF < 0.05 || lM < 0.05)              return("Concordant_single_arm")
    }
    return("Not_DEG")
  }
  if (opp_sign) {
    if (!is.na(lF) && !is.na(lM) && lF < 0.05 && lM < 0.05) return("Divergent")
    if ((!is.na(lF) && lF < 0.05) || (!is.na(lM) && lM < 0.05))
      return("Divergent_one_sided")
  }
  if (same_sign) {
    if (abs(bF) >= 2*abs(bM) && !is.na(lF) && lF < 0.05 && !is.na(lM)) {
      if (lM > 0.20) return("Female_biased")
      if (lM >= 0.05 && lM <= 0.20 && !is.na(pow) && pow < 0.5)
        return("Female_biased_M_underpowered")
    }
    if (abs(bM) >= 2*abs(bF) && !is.na(lM) && lM < 0.05 && !is.na(lF)) {
      if (lF > 0.20) return("Male_biased")
      if (lF >= 0.05 && lF <= 0.20) return("Male_biased_F_underpowered")
    }
    return("Sex_modifier")
  }
  return("Uncertain")
}
obs2[, class_v5_perm := mapply(classify_v5,
                                beta_F_strat, beta_M_strat,
                                q_emp_per_gene, lfsr_F_strat, lfsr_M_strat,
                                power_M_at_F, chr)]
obs2[chr == "chrY" & class_v5_perm == "Female_biased",
     class_v5_perm := "Male_biased"]

cat("\n== class_v5_perm distribution ==\n")
print(table(obs2$class_v5_perm, useNA = "ifany"))

# ---- v5 Tier-A retention check ----
if (file.exists(V5_CSV)) {
  v5 <- fread(V5_CSV, select = c("gene", "padj_int_5k", "evidence_tier_A_B_C"))
  v5_tierA <- v5[evidence_tier_A_B_C == "A", gene]
  perm_at_005 <- obs2[!is.na(q_emp_per_gene) & q_emp_per_gene < 0.05, gene]
  retain <- length(intersect(v5_tierA, perm_at_005)) / max(length(v5_tierA), 1)
  cat("v5 Tier-A retention at q_emp_per_gene<0.05:",
      round(retain * 100, 1), "% (", length(intersect(v5_tierA, perm_at_005)),
      "/", length(v5_tierA), ")\n")
}

# Counts
cat("\n== Perm-FDR counts ==\n")
cat("  q_emp_per_gene < 0.05:", sum(obs2$q_emp_per_gene < 0.05, na.rm = TRUE), "\n")
cat("  q_emp_per_gene < 0.10:", sum(obs2$q_emp_per_gene < 0.10, na.rm = TRUE), "\n")
cat("  q_emp_pooled   < 0.05:", sum(obs2$q_emp_pooled   < 0.05, na.rm = TRUE), "\n")
cat("  q_emp_pooled   < 0.10:", sum(obs2$q_emp_pooled   < 0.10, na.rm = TRUE), "\n")

out_cols <- c("gene", "t_int_obs", "p_dream_obs",
              "p_emp_per_gene", "q_emp_per_gene",
              "p_emp_pooled",   "q_emp_pooled",
              "n_perm_used", "n_exceed",
              "agree_dream_vs_perm", "class_v5_perm")
out <- obs2[, ..out_cols]
tmp <- paste0(OUT_CSV, ".tmp"); fwrite(out, tmp); file.rename(tmp, OUT_CSV)
cat("\nWrote:", OUT_CSV, "  rows:", nrow(out), "  cols:", ncol(out), "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
