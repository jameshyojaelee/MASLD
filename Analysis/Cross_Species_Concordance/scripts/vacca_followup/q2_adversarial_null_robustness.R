#!/usr/bin/env Rscript
# ============================================================================
# q2_adversarial_null_robustness.R   (Vacca benchmark — Q2 ADVERSARIAL)
#
# Confirm the threshold-free concordance findings are NOT null artifacts of
# set size / ranking structure, via explicit LABEL-SHUFFLE permutation nulls
# (1000x each), independent of fgsea's internal gene-set permutation.
#
#   (A) q2a GSEA null  — observed: Vacca DSEA-951 signature is enriched at the
#       top of our continuous concordance rankings (translatability_score
#       NES 1.445; mean_h_lfc NES 2.42; q2_threshold_free_fgsea.csv).
#       NULL: draw a random gene set of the SAME size (= #Vacca-in-universe)
#       1000x and recompute the Kolmogorov–Smirnov-style running enrichment
#       score (ES) for each ranking. Empirical p = fraction of null ES >= obs ES.
#       (i.e. shuffle which genes are "Vacca members" → the membership-shuffle.)
#
#   (B) q2c correlation null — observed: our cutoff-free per-gene concordance
#       metrics correlate with Vacca's continuous progression logFC
#       (translatability_score vs |vacca_prog_lfc| Spearman 0.220, n=11178;
#        hm_agreement vs signed vacca_prog_lfc Spearman 0.211, n=11223;
#        q2_continuous_concordance_vs_vacca_prog.csv).
#       NULL: permute the gene<->vacca_prog_lfc pairing 1000x and recompute
#       Spearman. Empirical p = fraction of |null rho| >= |obs rho|.
#
# Reuses q2_pergene_merged.csv (exact observed-correlation inputs) for (B) and
# rebuilds the atlas rankings + Vacca-951 signature for (A), matching the
# original scripts' construction exactly.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
set.seed(42)
N_PERM <- 1000L
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# (A) q2a GSEA membership-shuffle null
# ============================================================================
# --- Vacca 951-gene DSEA signature (identical construction to q2a script) ---
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet = "Table S4"))
dsea_col <- grep("Gene_used_in_DSEA", names(s4), value = TRUE)[1]
sym_col  <- grep("^GeneSymbol$", names(s4), value = TRUE)[1]
if (is.na(sym_col)) sym_col <- grep("GeneSymbol", names(s4), value = TRUE)[1]
s4[, dsea := suppressWarnings(as.numeric(get(dsea_col)))]
vacca_sig <- unique(toupper(trimws(s4[dsea == 1][[sym_col]])))
vacca_sig <- vacca_sig[!is.na(vacca_sig) & vacca_sig != ""]

# --- atlas rankings (identical collapse/ordering to q2a script) ---
a <- fread(file.path(CS, "concordance_atlas_unified.csv"))
a[, sym := toupper(trimws(human_symbol))]
a <- a[!is.na(sym) & sym != ""]
setorder(a, -n_concordant, -translatability_score)
a <- a[!duplicated(sym)]
a[, signed_nconc := sign(mean_h_lfc) * n_concordant]
in_universe <- intersect(vacca_sig, a$sym)
set_size <- length(in_universe)
cat(sprintf("(A) atlas genes=%d ; Vacca-951 in universe (set_size)=%d\n",
            nrow(a), set_size))

# weighted (p=1) running-sum enrichment score: classic GSEA ES for a ranked
# stats vector and a membership index set. Returns the signed max-deviation ES.
gsea_ES <- function(stats_sorted, hit_idx) {
  N <- length(stats_sorted)
  hit <- logical(N); hit[hit_idx] <- TRUE
  w <- abs(stats_sorted)
  Nh <- sum(w[hit])
  if (Nh == 0) return(0)
  step <- numeric(N)
  step[hit]  <-  w[hit] / Nh
  step[!hit] <- -1 / (N - length(hit_idx))
  rs <- cumsum(step)
  rs[which.max(abs(rs))]
}

make_rank <- function(dt, statcol) {
  v <- dt[[statcol]]
  jit <- (seq_len(nrow(dt)) / nrow(dt)) * 1e-6
  setNames(v + jit, dt$sym)
}

rankings <- c("translatability_score", "mean_h_lfc", "signed_nconc", "n_concordant")
resA <- rbindlist(lapply(rankings, function(rk) {
  stats <- sort(make_rank(a, rk), decreasing = TRUE)
  nm    <- names(stats)
  obs_idx <- which(nm %in% in_universe)
  obs_ES  <- gsea_ES(stats, obs_idx)
  # membership-shuffle null: random gene sets of the same size
  Nall <- length(stats)
  null_ES <- numeric(N_PERM)
  for (p in seq_len(N_PERM)) {
    idx <- sample.int(Nall, set_size)
    null_ES[p] <- gsea_ES(stats, idx)
  }
  # one-sided test in the OBSERVED direction (these rankings give obs_ES>0)
  if (obs_ES >= 0) {
    emp_p <- (1 + sum(null_ES >= obs_ES)) / (N_PERM + 1)
  } else {
    emp_p <- (1 + sum(null_ES <= obs_ES)) / (N_PERM + 1)
  }
  z <- (obs_ES - mean(null_ES)) / sd(null_ES)
  data.table(test = "q2a_gsea_membershuffle", ranking = rk,
             obs_ES = round(obs_ES, 4),
             null_ES_mean = round(mean(null_ES), 4),
             null_ES_sd = round(sd(null_ES), 4),
             null_ES_max = round(max(null_ES), 4),
             z = round(z, 3),
             emp_p = emp_p, n_perm = N_PERM, set_size = set_size)
}))
cat("\n=== (A) q2a GSEA membership-shuffle null (1000x random same-size sets) ===\n")
print(resA)

# ============================================================================
# (B) q2c continuous-correlation label-shuffle null
# ============================================================================
prim <- fread(file.path(OUT, "q2_pergene_merged.csv"))
# Recompute the two PRIMARY observed cutoff-free concordances, exactly as q2c:
#   (1) translatability_score vs abs(vacca_prog_lfc)
#   (2) hm_agreement vs vacca_prog_lfc (signed)
# plus the two extra continuous metrics reported, for completeness.
sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))

perm_spearman_p <- function(x, target, two_sided = TRUE) {
  ok <- is.finite(x) & is.finite(target)
  x <- x[ok]; target <- target[ok]; n <- length(x)
  obs <- sp(x, target)
  null <- numeric(N_PERM)
  for (p in seq_len(N_PERM)) null[p] <- sp(x, target[sample.int(n)])
  if (two_sided) {
    emp_p <- (1 + sum(abs(null) >= abs(obs))) / (N_PERM + 1)
  } else {
    emp_p <- (1 + sum(null >= obs)) / (N_PERM + 1)
  }
  z <- (obs - mean(null)) / sd(null)
  list(obs = obs, null_mean = mean(null), null_sd = sd(null),
       null_abs_max = max(abs(null)), z = z, emp_p = emp_p, n = n)
}

corr_specs <- list(
  list(metric = "translatability_score", target = "abs_vacca_prog_lfc",
       x = prim$translatability_score, t = abs(prim$vacca_prog_lfc)),
  list(metric = "hm_agreement", target = "vacca_prog_lfc_signed",
       x = prim$hm_agreement, t = prim$vacca_prog_lfc),
  list(metric = "translatability_score", target = "vacca_prog_lfc_signed",
       x = prim$translatability_score, t = prim$vacca_prog_lfc),
  list(metric = "hm_agreement", target = "abs_vacca_prog_lfc",
       x = prim$hm_agreement, t = abs(prim$vacca_prog_lfc)),
  list(metric = "atlas_mean_h_lfc", target = "vacca_prog_lfc_signed",
       x = prim$mean_h_lfc, t = prim$vacca_prog_lfc),
  list(metric = "human_ss_only", target = "vacca_prog_lfc_signed",
       x = prim$human_ss, t = prim$vacca_prog_lfc)
)

resB <- rbindlist(lapply(corr_specs, function(s) {
  r <- perm_spearman_p(s$x, s$t, two_sided = TRUE)
  data.table(test = "q2c_corr_labelshuffle", metric = s$metric, target = s$target,
             obs_spearman = round(r$obs, 4),
             null_mean = round(r$null_mean, 4),
             null_sd = round(r$null_sd, 4),
             null_abs_max = round(r$null_abs_max, 4),
             z = round(r$z, 2),
             emp_p = r$emp_p, n = r$n, n_perm = N_PERM)
}))
cat("\n=== (B) q2c continuous-correlation label-shuffle null (1000x permuted pairing) ===\n")
print(resB)

# ============================================================================
# Combine + write
# ============================================================================
combined <- rbindlist(list(
  resA[, .(test, name = ranking, observed = obs_ES, null_center = null_ES_mean,
           null_sd = null_ES_sd, null_extreme = null_ES_max, z, emp_p, n = set_size)],
  resB[, .(test, name = paste0(metric, "__", target), observed = obs_spearman,
           null_center = null_mean, null_sd = null_sd, null_extreme = null_abs_max,
           z, emp_p, n)]
), use.names = TRUE)
fwrite(resA, file.path(OUT, "q2_adversarial_null_gsea.csv"))
fwrite(resB, file.path(OUT, "q2_adversarial_null_corr.csv"))
fwrite(combined, file.path(OUT, "q2_adversarial_null_robustness.csv"))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("VERDICT (empirical p over", N_PERM, "label/membership shuffles):\n")
cat(sprintf("  q2a translatability_score GSEA : obs_ES=%.3f  emp_p=%.4g  (null max ES=%.3f)\n",
            resA[ranking == "translatability_score", obs_ES],
            resA[ranking == "translatability_score", emp_p],
            resA[ranking == "translatability_score", null_ES_max]))
cat(sprintf("  q2a mean_h_lfc GSEA            : obs_ES=%.3f  emp_p=%.4g\n",
            resA[ranking == "mean_h_lfc", obs_ES],
            resA[ranking == "mean_h_lfc", emp_p]))
cat(sprintf("  q2c translatability vs |prog|  : rho=%.3f  emp_p=%.4g  (null |rho| max=%.3f)\n",
            resB[metric == "translatability_score" & target == "abs_vacca_prog_lfc", obs_spearman],
            resB[metric == "translatability_score" & target == "abs_vacca_prog_lfc", emp_p],
            resB[metric == "translatability_score" & target == "abs_vacca_prog_lfc", null_abs_max]))
cat(sprintf("  q2c hm_agreement vs prog(signed): rho=%.3f  emp_p=%.4g\n",
            resB[metric == "hm_agreement" & target == "vacca_prog_lfc_signed", obs_spearman],
            resB[metric == "hm_agreement" & target == "vacca_prog_lfc_signed", emp_p]))
cat(strrep("=", 78), "\n")
