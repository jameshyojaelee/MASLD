#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# T1 GATE-B2 / GATE-N2 — NB spike-in positive control + mean-variance decoupling.
#
# The whole T1 claim ("variance-only genes invisible to logFC") is valid ONLY if
# the DV statistic (a) HAS POWER on true variance shifts and (b) does NOT flag
# genes whose only change is a MEAN shift (NB mean-variance coupling could leak).
# Synthetic NB counts with 4 gene classes — null / mean-only / var-only / both —
# must yield: high DV power on {var,both}; and after the decoupling gate
# (residualise dv_t on |mean logFC|), mean-only genes fall back to the null.
# Lightweight (toy NB data) — safe off-node.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/heterogeneity_program/lib/dv_stat.R"))

set.seed(7)
nb <- 2L; per_grp <- 50L                       # 2 datasets × (50 ctrl + 50 dis) = 200
block <- factor(rep(paste0("ds", seq_len(nb)), each = 2L * per_grp))
group <- factor(rep(rep(c("ctrl", "dis"), each = per_grp), times = nb), levels = c("ctrl","dis"))
N <- length(group); dis <- group == "dis"

G <- 600L
cls <- rep(c("null","mean","var","both"), times = c(300L, 100L, 100L, 100L))
mu0  <- 2^runif(G, 4, 11)                       # baseline mean counts
phi0 <- runif(G, 0.1, 0.5)                      # baseline NB dispersion
fc   <- ifelse(cls %in% c("mean","both"), 2^(sample(c(-1,1),G,TRUE) * runif(G,0.8,1.6)), 1)  # mean shift
kv   <- ifelse(cls %in% c("var","both"),  runif(G, 2.5, 4.0), 1)                              # dispersion ×
libf <- runif(N, 0.6, 1.6)                       # per-sample library factor

counts <- matrix(0L, G, N)
for (j in seq_len(N)) {
  mu  <- mu0 * (if (dis[j]) fc else 1) * libf[j]
  phi <- phi0 * (if (dis[j]) kv else 1)
  counts[, j] <- rnbinom(G, size = 1/phi, mu = mu)
}
rownames(counts) <- sprintf("g%04d", seq_len(G))
dge <- DGEList(counts); dge <- calcNormFactors(dge)

# observed DV (full sample), fixed gene universe
keep <- filterByExpr(dge, group = group)
gs <- rownames(dge)[keep]
res <- dv_screen(dge, group, block, gene_set = gs)
res <- as.data.table(res); res[, cls := cls[match(gene, rownames(counts))]]

# "mean logFC" proxy = observed disease/control log-mean (the mean-shift axis)
logcpm <- edgeR::cpm(dge, log = TRUE)
res[, abs_logfc := abs(rowMeans(logcpm[gene, dis]) - rowMeans(logcpm[gene, !dis]))]

# Decoupling is INTRINSIC to limma-DV on log-CPM (no global residualisation —
# that subtracts real signal from the "both" genes). It is demonstrated by:
#  (i) mean-only genes not flagged, (ii) corr(dv_t,|logFC|)≈0 AMONG genes with no
#  true variance change, (iii) the variance-only headline restricted to low-|logFC|.
q <- function(p) p.adjust(p, "BH")
res[, dv_q := q(dv_p)]
nonvar <- res$cls %in% c("null", "mean")          # NO true variance change
corr_all    <- res[, cor(dv_t, abs_logfc)]        # inflated by legitimate "both" genes
corr_nonvar <- res[nonvar, cor(dv_t, abs_logfc)]  # the real artifact check

pw <- function(cc) res[cls == cc, mean(dv_q < 0.10)]
cat("── DV power / false-positive by gene class ──\n")
cat(sprintf("  var-only  power (q<0.10)        = %.2f   (target high)\n", pw("var")))
cat(sprintf("  both      power (q<0.10)        = %.2f   (target high)\n", pw("both")))
cat(sprintf("  mean-only FLAGGED (q<0.10)      = %.2f   (target ≈ null rate)\n", pw("mean")))
cat(sprintf("  null      FLAGGED (q<0.10)      = %.2f   (target ≈0)\n", pw("null")))
cat(sprintf("  corr(dv_t,|logFC|) ALL          = %+.3f   (inflated by real 'both' genes)\n", corr_all))
cat(sprintf("  corr(dv_t,|logFC|) NON-variance = %+.3f   (target ≈0 — the artifact check)\n", corr_nonvar))

pass <- pw("var") > 0.6 && pw("both") > 0.6 &&
        pw("mean") < 0.15 && pw("null") < 0.15 && abs(corr_nonvar) < 0.15
cat("\n=== ", if (pass) "PASS" else "FAIL",
    " (DV powered on variance; mean shifts do NOT leak into the no-variance set) ===\n", sep = "")
if (!pass) quit(status = 1)
