#!/usr/bin/env Rscript
# BF2: continuum position at which each molecular system reaches half its total difference (PRESPEC_BF2.md).
suppressPackageStartupMessages({library(mgcv); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]; B <- as.integer(ifelse(length(args) > 1, args[2], 1000))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CG0 <- file.path(ROOT, "Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z")
META <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/molecular_systems_inference",
                  "hac-system-inference-corrected774-20261009T171017Z/tables/system_continuum_meta_analysis.tsv")
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
set.seed(20261023)
RNGkind("L'Ecuyer-CMRG")

u <- read.delim(file.path(CG0, "units.tsv"))
meta <- read.delim(META)
meta <- meta[meta$aggregation_method == "collection_balanced" & meta$axis_id == "signature_pc1", ]
eligible <- sort(meta$community_id[meta$random_q_value < 0.05])
systems <- sprintf("sys_%02d", 1:43)
elig_cols <- sub("system_", "sys_", eligible)
cat("eligible systems:", length(eligible), "\n")

within_stage <- function(d) {
  d <- d[!is.na(d$fibrosis_stage), ]
  d[!(d$dataset == "GSE213621" & d$fibrosis_stage == 0), ]
}
sets <- list(
  GSE162694 = u[u$dataset == "GSE162694", ],
  GSE213621 = u[u$dataset == "GSE213621", ],
  GSE162694_within_stage = within_stage(u[u$dataset == "GSE162694", ]),
  GSE213621_within_stage = within_stage(u[u$dataset == "GSE213621", ]),
  GSE130970 = u[u$dataset == "GSE130970", ],
  GSE135251 = u[u$dataset == "GSE135251", ])
stopifnot(nrow(sets$GSE162694) == 140, nrow(sets$GSE213621) == 326,
          nrow(sets$GSE162694_within_stage) == 109, nrow(sets$GSE213621_within_stage) == 258)

PCT <- seq(2.5, 97.5, by = 0.5)

curve_stats <- function(d, col) {
  d$y <- d[[col]]
  d <- d[!is.na(d$y) & !is.na(d$continuum_z), ]
  d$inferred_sex <- factor(d$inferred_sex)
  two_sex <- length(unique(d$inferred_sex)) > 1
  f <- if (two_sex) y ~ s(continuum_z, k = 6, bs = "cr") + inferred_sex else y ~ s(continuum_z, k = 6, bs = "cr")
  fit <- gam(f, data = d, method = "REML")
  x <- quantile(d$continuum_z, PCT / 100, type = 7, names = FALSE)
  nd <- data.frame(continuum_z = x, inferred_sex = factor(levels(d$inferred_sex)[1], levels = levels(d$inferred_sex)))
  fx <- predict(fit, newdata = nd, type = "terms")[, "s(continuum_z)"]
  delta <- fx[length(fx)] - fx[1]; s <- sign(delta)
  prog <- s * (fx - fx[1])
  at <- function(fr) PCT[min(which(prog >= fr * abs(delta)))]
  lin <- fx[1] + (fx[length(fx)] - fx[1]) * (x - x[1]) / (x[length(x)] - x[1])
  dev <- s * (fx - lin)  # positive: curve is ahead of the straight line (concave toward its end level)
  list(stats = unlist(lapply(list(delta = delta, delta_sd = delta / sd(d$y), m = at(0.5), q25 = at(0.25), q75 = at(0.75),
                 max_signed_dev = dev[which.max(abs(dev))] / abs(delta), edf = sum(fit$edf[grepl("s\\(", names(fit$edf))]),
                 p_smooth = summary(fit)$s.table[1, "p-value"]), unname)),
       curve = fx)
}

point <- list(); curves <- list()
for (nm in names(sets)) {
  for (col in systems) {
    r <- curve_stats(sets[[nm]], col)
    point[[length(point) + 1]] <- data.frame(set = nm, system = sub("sys_", "system_", col), t(r$stats), n = nrow(sets[[nm]]))
    if (nm %in% c("GSE162694", "GSE213621")) curves[[length(curves) + 1]] <- data.frame(set = nm, system = sub("sys_", "system_", col), pct = PCT, f = r$curve)
  }
}
point <- do.call(rbind, point); point$eligible <- point$system %in% eligible
write.table(point, file.path(out, "midpoints_point.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, curves), file.path(out, "curves_primary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# bootstrap: primary cohorts, eligible systems
boot <- list()
for (nm in c("GSE162694", "GSE213621")) {
  d <- sets[[nm]]
  res <- mclapply(seq_len(B), function(b) {
    db <- d[sample.int(nrow(d), replace = TRUE), ]
    sapply(elig_cols, function(col) tryCatch(curve_stats(db, col)$stats["m"], error = function(e) NA_real_))
  }, mc.cores = cores, mc.set.seed = TRUE)
  M <- do.call(rbind, res); colnames(M) <- eligible
  boot[[nm]] <- M
  write.table(data.frame(draw = seq_len(B), M, check.names = FALSE), file.path(out, paste0("boot_m_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# primary test
pm <- function(nm) { x <- point[point$set == nm & point$eligible, ]; setNames(x$m, x$system)[eligible] }
m1 <- pm("GSE162694"); m2 <- pm("GSE213621")
rho <- cor(m1, m2, method = "spearman")
perm <- replicate(10000, cor(m1, sample(m2), method = "spearman"))
p_perm <- (1 + sum(abs(perm) >= abs(rho))) / 10001
rb <- sapply(seq_len(B), function(b) cor(boot$GSE162694[b, ], boot$GSE213621[b, ], method = "spearman", use = "complete.obs"))
ci <- quantile(rb, c(0.025, 0.975), na.rm = TRUE)
bci <- do.call(rbind, lapply(names(boot), function(nm) data.frame(set = nm, system = eligible,
  m_lo = apply(boot[[nm]], 2, quantile, 0.025, na.rm = TRUE), m_hi = apply(boot[[nm]], 2, quantile, 0.975, na.rm = TRUE),
  frac_na = colMeans(is.na(boot[[nm]])))))
write.table(bci, file.path(out, "midpoints_boot_ci.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sink(file.path(out, "primary_test.txt"))
cat(sprintf("n_eligible\t%d\nrho\t%.4f\np_perm\t%.5f\nboot_ci_lo\t%.4f\nboot_ci_hi\t%.4f\ndecision_order_replicates\t%s\n",
            length(eligible), rho, p_perm, ci[1], ci[2], (p_perm < 0.05) & (ci[1] > 0)))
sink()
cat(readLines(file.path(out, "primary_test.txt")), sep = "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
