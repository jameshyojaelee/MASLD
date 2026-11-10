# RV5: per-gene heterogeneity of the 1,854 bulk calls across five cohorts. Follows PRESPEC_RV5.md.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
COH <- file.path(ROOT, "figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv")
DISC <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
wr <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")

disc <- fread(DISC)[, .(gene = gene, symbol, d_logFC = logFC, call = padj < 0.05 & abs(logFC) > 0.5)]
stopifnot(sum(disc$call) == 1854)
r <- fread(COH, select = c("gene_id_versioned", "dataset", "logFC", "t"))
r[, se := logFC / t]
r <- r[is.finite(logFC) & is.finite(se) & se > 0]
cohorts <- sort(unique(r$dataset)); stopifnot(length(cohorts) == 5)
B <- as.matrix(dcast(r, gene_id_versioned ~ dataset, value.var = "logFC"), rownames = "gene_id_versioned")[, cohorts]
S <- as.matrix(dcast(r, gene_id_versioned ~ dataset, value.var = "se"), rownames = "gene_id_versioned")[, cohorts]

re_meta <- function(b, s) {
  w <- 1 / s^2; k <- rowSums(!is.na(b))
  W <- rowSums(w, na.rm = TRUE); fe <- rowSums(w * b, na.rm = TRUE) / W
  Q <- rowSums(w * (b - fe)^2, na.rm = TRUE); df <- k - 1
  C <- W - rowSums(w^2, na.rm = TRUE) / W
  tau2 <- pmax(0, (Q - df) / C)
  wr_ <- 1 / (s^2 + tau2); Wr <- rowSums(wr_, na.rm = TRUE)
  mu <- rowSums(wr_ * b, na.rm = TRUE) / Wr; se <- sqrt(1 / Wr)
  data.table(k = k, fe = fe, Q = Q, tau = sqrt(tau2), mu = mu, se_re = se, p_z = 2 * pnorm(-abs(mu / se)),
             I2 = ifelse(Q > 0, pmax(0, (Q - df) / Q), 0), p_Q = pchisq(Q, df, lower.tail = FALSE),
             pi_lo = mu - qt(0.975, pmax(k - 2, 1)) * sqrt(tau2 + se^2), pi_hi = mu + qt(0.975, pmax(k - 2, 1)) * sqrt(tau2 + se^2))
}
bh_k3 <- function(p, k) { q <- rep(NA_real_, length(p)); i <- k >= 3; q[i] <- p.adjust(p[i], "BH"); q }

m <- re_meta(B, S); m[, gene := rownames(B)]
m[, q_z := bh_k3(p_z, k)]; m[, q_Q := bh_k3(p_Q, k)]
loo <- lapply(cohorts, function(c) {
  keep <- setdiff(cohorts, c)
  x <- re_meta(B[, keep], S[, keep])
  data.table(gene = rownames(B), mu = x$mu, q = bh_k3(x$p_z, x$k))
})
names(loo) <- cohorts

h <- merge(disc[call == TRUE], m, by = "gene")
stopifnot(nrow(h) == 1854)
sg <- sign(h$d_logFC)
Bc <- B[h$gene, , drop = FALSE]
h[, n_opposite := rowSums(sign(Bc) == -sg, na.rm = TRUE)]
h[, n_same := rowSums(sign(Bc) == sg, na.rm = TRUE)]
for (c in cohorts) {
  l <- loo[[c]][match(h$gene, gene)]
  h[, paste0("loo_mu_", c) := l$mu]; h[, paste0("loo_q_", c) := l$q]
}
loo_fail <- sapply(cohorts, function(c) {
  mu <- h[[paste0("loo_mu_", c)]]; q <- h[[paste0("loo_q_", c)]]
  is.na(q) | q >= 0.05 | sign(mu) != sg
})
h[, loo_fail_cohorts := apply(loo_fail, 1, function(x) paste(cohorts[x], collapse = ","))]
re_ok <- !is.na(h$q_z) & h$q_z < 0.05 & sign(h$mu) == sg
dir_ok <- h$k >= 4 & h$n_opposite <= 1
h[, class := fifelse(!dir_ok, "cohort_discordant",
              fifelse(!re_ok, "direction_consistent_size_heterogeneous",
              fifelse(loo_fail_cohorts == "", "consistent", "single_cohort_dependent")))]
h[, pi_excludes_0 := pi_lo > 0 | pi_hi < 0]
h[, direction := fifelse(d_logFC > 0, "up", "down")]
wr(h, "call_heterogeneity.tsv")

counts <- h[, .N, by = .(class, direction)][order(class, direction)]
wr(counts, "class_counts.tsv")
sink(file.path(OUT, "summary.txt"))
cat("Cohorts:", paste(cohorts, collapse = ", "), "\n\nClass counts\n"); print(h[, .N, by = class][order(-N)]); print(dcast(counts, class ~ direction, value.var = "N", fill = 0))
cat("\nSingle-cohort dependent: cohort(s) whose omission loses the call\n"); print(h[class == "single_cohort_dependent", .N, by = loo_fail_cohorts][order(-N)])
cat("\nMedian tau and I2 by class\n"); print(h[, .(median_tau = median(tau), median_I2 = median(I2), Q_q05 = sum(q_Q < 0.05, na.rm = TRUE)), by = class])
cat(sprintf("\nPrediction interval excludes 0: %d of 1854\n", sum(h$pi_excludes_0, na.rm = TRUE)))
cat(sprintf("Cochran Q BH q < 0.05: %d; RE z BH q < 0.05 same sign: %d\n", sum(h$q_Q < 0.05, na.rm = TRUE), sum(re_ok)))
cat("\nOpposite-sign cohort count\n"); print(table(h$n_opposite, h$k))
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done\n")
