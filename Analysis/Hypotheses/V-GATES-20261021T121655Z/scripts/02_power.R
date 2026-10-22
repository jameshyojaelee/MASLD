# WP0d: minimum detectable effect (SD units) and power for each planned design.
# Two-group designs: noncentral t, df = n1 + n2 - 2 - k (k extra covariates).
# Paired / sign-flip designs: one-sample noncentral t, df = n - 1, plus the
# smallest attainable exact one-sided sign-flip p (1 / 2^n).
# Usage: Rscript 02_power.R <run_dir>
args <- commandArgs(trailingOnly = TRUE)
out_dir <- file.path(args[1], "design")
alpha <- 0.05

power_two <- function(d, n1, n2, k) {
  df <- n1 + n2 - 2 - k
  ncp <- d / sqrt(1 / n1 + 1 / n2)
  tc <- qt(1 - alpha / 2, df)
  pt(tc, df, ncp, lower.tail = FALSE) + pt(-tc, df, ncp)
}
power_one <- function(d, n) {
  df <- n - 1
  ncp <- d * sqrt(n)
  tc <- qt(1 - alpha / 2, df)
  pt(tc, df, ncp, lower.tail = FALSE) + pt(-tc, df, ncp)
}
mde <- function(f) uniroot(function(d) f(d) - 0.80, c(0.01, 20))$root

designs <- list(
  list(id = "V5 GSE298719 disease vs healthy", type = "two", n1 = 6, n2 = 7, k = 2),
  list(id = "V5 GSE298719, no covariates", type = "two", n1 = 6, n2 = 7, k = 0),
  list(id = "V5 within-plate subset (E35250325_1 + E36240625_1)", type = "two", n1 = 5, n2 = 5, k = 1),
  list(id = "V5 without F4 donors", type = "two", n1 = 6, n2 = 5, k = 2),
  list(id = "V8 GSE344087 fibrosis yes vs no", type = "two", n1 = 16, n2 = 3, k = 0),
  list(id = "V7d SCP3359 MASH cirrhosis vs control", type = "two", n1 = 9, n2 = 9, k = 0),
  list(id = "V6 GSE255312 within-ROI, all ROIs", type = "one", n = 24),
  list(id = "V7d SCP3359 within-donor gradient, cirrhotic", type = "one", n = 9),
  list(id = "V7c GSE292268 within-biopsy", type = "one", n = 4),
  list(id = "V7b Fukumoto within-donor", type = "one", n = 2)
)
ref_d <- c(discovery_projection = 0.61, v3_c2 = 1.08, large = 1.5)

rows <- lapply(designs, function(x) {
  if (x$type == "two") {
    f <- function(d) power_two(d, x$n1, x$n2, x$k)
    n_lab <- sprintf("%d v %d", x$n1, x$n2); minp <- NA
  } else {
    f <- function(d) power_one(d, x$n)
    n_lab <- sprintf("%d paired", x$n); minp <- 1 / 2^x$n
  }
  m <- if (x$type == "one" && x$n < 3) NA else mde(f)
  pw <- if (x$type == "one" && x$n < 3) rep(NA, 3) else vapply(ref_d, f, numeric(1))
  data.frame(design = x$id, n = n_lab, mde_sd_80pct = round(m, 2),
             power_d0.61 = round(pw[1], 2), power_d1.08 = round(pw[2], 2),
             power_d1.50 = round(pw[3], 2), min_exact_signflip_p_one_sided = signif(minp, 3))
})
tab <- do.call(rbind, rows)
write.table(tab, file.path(out_dir, "mde_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(tab, row.names = FALSE)

# Pooled-test type I error conditioned on V3's observed z (1.83), sqrt(n) weights.
z1 <- 1.83; w1 <- sqrt(17); w2 <- sqrt(13)
z2_needed <- (qnorm(0.975) * sqrt(w1^2 + w2^2) - w1 * z1) / w2
cat(sprintf("\nPooled V3+V5 test: V5 z needed %.3f; conditional type I error %.3f\n",
            z2_needed, pnorm(z2_needed, lower.tail = FALSE)))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_power.txt"))
