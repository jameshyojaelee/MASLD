# U8a: independent within-dataset permutation null for H59 (Fig 4C).
# Usage: Rscript --vanilla 02_null.R <seed> <B> <outdir> <ncores>
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1]); B <- as.integer(args[2]); out <- args[3]; ncores <- as.integer(args[4])
here <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive"
source(file.path(here, "u8a_lib.R"))
suppressPackageStartupMessages(library(parallel))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()

u62 <- load_donor_scores()
u60 <- u62[dataset != "GSE189600"]
fit <- u60[dataset %in% FIT_DS]
progs <- sort(unique(fit$program))
donors62 <- sort(unique(u62$donor), method = "radix")
stopifnot(length(donors62) == 62L, length(progs) == 117L)
K <- make_keys(seed, B, donors62)
say("seed %d B %d keys %d x %d; K[1,1:3] = %s", seed, B, nrow(K), ncol(K), paste(signif(K[1, 1:3], 8), collapse = " "))

# closed-form equal-variance OLS t with dataset intercepts
ols_t <- function(y, s, d) {
  yc <- y - ave(y, d); sc <- s - ave(s, d)
  sxx <- sum(sc^2); b <- sum(sc * yc) / sxx
  df <- length(y) - length(unique(d)) - 1
  rss <- sum((yc - b * sc)^2)
  b / sqrt(rss / df / sxx)
}

obs <- rbindlist(lapply(progs, function(p) {
  d <- fit[program == p]; g <- gls_fit(d$score, d$stage, d$dataset)
  e62 <- u62[program == p]; e60 <- u60[program == p]
  data.table(program = p, t = g$t, beta = g$beta, se = g$se,
             t62 = ols_t(e62$score, e62$stage, e62$dataset), t60 = ols_t(e60$score, e60$stage, e60$dataset))
}))

one_prog <- function(p) {
  d <- fit[program == p]; e62 <- u62[program == p]; e60 <- u60[program == p]
  tg <- t62 <- t60 <- numeric(B)
  for (k in seq_len(B)) {
    kr <- K[k, ]
    sp <- permute_stage(d$donor, d$stage, d$dataset, kr)
    tg[k] <- tryCatch(gls_fit(d$score, sp, d$dataset)$t, error = function(e) NA_real_)
    t62[k] <- ols_t(e62$score, permute_stage(e62$donor, e62$stage, e62$dataset, kr), e62$dataset)
    t60[k] <- ols_t(e60$score, permute_stage(e60$donor, e60$stage, e60$dataset, kr), e60$dataset)
  }
  list(tg = tg, t62 = t62, t60 = t60)
}
nl <- mclapply(progs, one_prog, mc.cores = ncores, mc.preschedule = FALSE)
names(nl) <- progs
TG <- sapply(nl, `[[`, "tg"); T62 <- sapply(nl, `[[`, "t62"); T60 <- sapply(nl, `[[`, "t60")
say("null done in %.1f min; failed GLS draws %d", as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(!is.finite(TG)))
fwrite(as.data.table(TG), file.path(out, "null_gls_t.tsv.gz"), sep = "\t")
fwrite(as.data.table(T62), file.path(out, "null_ols62_t.tsv.gz"), sep = "\t")

reg <- fread(REGISTRY_FILE); reg[, program := paste0(cell_type, "_", module)]
obs <- merge(obs, reg[, .(program, module_name, stability_median)], by = "program")
obs[, perm_p := sapply(seq_len(.N), function(i) perm_p(t[i], TG[, program[i]]))]
obs[, perm_q := p.adjust(perm_p, "BH", n = 117)]
obs[, perm_p62 := sapply(seq_len(.N), function(i) perm_p(t62[i], T62[, program[i]]))]
obs[, perm_q62 := p.adjust(perm_p62, "BH", n = 117)]
obs[, perm_p60 := sapply(seq_len(.N), function(i) perm_p(t60[i], T60[, program[i]]))]
obs[, perm_q60 := p.adjust(perm_p60, "BH", n = 117)]
obs[, crit := sapply(program, function(p) quantile(abs(TG[, p]), 1 - 0.05 / 117, type = 7, na.rm = TRUE))]
obs[, selected := perm_q < 0.05 & stability_median >= 0.5]
# maxT family-wise p
mx <- apply(abs(TG), 1, max, na.rm = TRUE)
obs[, maxT_p := sapply(t, function(x) (1 + sum(mx >= abs(x))) / (1 + length(mx)))]
setorder(obs, perm_p)
fwrite(obs, file.path(out, "perm_results_117.tsv"), sep = "\t")
say("selected %d; q<0.05 %d; smallest q %.6f (%s); smallest maxT p %.5f", sum(obs$selected), sum(obs$perm_q < 0.05),
    min(obs$perm_q), obs$program[which.min(obs$perm_q)], min(obs$maxT_p))
print(obs[1:8, .(program, module_name, t, perm_p, perm_q, crit, perm_q62)])
print(obs[program %in% c("hepatocytes_20", "hepatocytes_8"), .(program, t, perm_p, perm_q, perm_q62, perm_q60, crit)])

# NC1: first 2,000 draws as pseudo-observed; BH over 117
nk <- min(2000L, B)
A <- abs(TG)
calls <- integer(nk)
srt <- lapply(seq_len(ncol(A)), function(j) sort(A[, j]))
for (k in seq_len(nk)) {
  pk <- sapply(seq_len(ncol(A)), function(j) {
    x <- A[k, j]; ge <- length(srt[[j]]) - findInterval(x, srt[[j]], left.open = TRUE)
    (1 + ge - 1) / B
  })
  calls[k] <- sum(p.adjust(pk, "BH", n = 117) < 0.05)
}
say("NC1: fraction of first %d draws with >=1 BH call: %.4f; max calls %d", nk, mean(calls > 0), max(calls))

# PC4 detection for the two named programs, rows 1..1000 of K
pc4 <- rbindlist(lapply(list(c("hepatocytes_20", 0.402086), c("hepatocytes_8", 0.351525)), function(z) {
  p <- z[1]; d <- fit[program == p]; cr <- obs[program == p]$crit
  rbindlist(lapply(c(0.2, 0.3, as.numeric(z[2])), function(delta) {
    det <- sapply(seq_len(min(1000L, B)), function(k) {
      sp <- permute_stage(d$donor, d$stage, d$dataset, K[k, ])
      abs(gls_fit(d$score + delta * sp, sp, d$dataset)$t) > cr
    })
    data.table(program = p, delta = delta, crit = cr, detection = mean(det))
  }))
}))
print(pc4)
fwrite(pc4, file.path(out, "pc4.tsv"), sep = "\t")

# PC3: own known-truth simulation (own seed), GLS 95% CI coverage
set.seed(seed + 1L)
des <- data.table(dataset = rep(FIT_DS, c(37, 18)), stage = c(rep(0:2, c(4, 7, 26)), rep(0:2, c(5, 4, 9))))
sdv <- ifelse(des$dataset == FIT_DS[1], 0.33, 0.95)
pc3 <- rbindlist(lapply(c(0, 0.30), function(b) {
  r <- t(sapply(1:2000, function(i) {
    y <- b * des$stage + rnorm(nrow(des), 0, sdv)
    g <- gls_fit(y, des$stage, des$dataset)
    h <- qt(0.975, 52) * g$se
    c(cover = (g$beta - h <= b) & (b <= g$beta + h), rej = abs(g$t) > qt(0.975, 52))
  }))
  data.table(true_beta = b, coverage = mean(r[, 1]), rejection = mean(r[, 2]))
}))
print(pc3)
fwrite(pc3, file.path(out, "pc3.tsv"), sep = "\t")

writeLines(c(sprintf("seed\t%d", seed), sprintf("B\t%d", B), "family\t117 (BH n=117)",
             sprintf("runtime_min\t%.1f", as.numeric(difftime(Sys.time(), t0, units = "mins")))), file.path(out, "run_params.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
system2("sha256sum", c(SCORE_FILE, META_FILE, REGISTRY_FILE), stdout = file.path(out, "inputs.sha256"))
