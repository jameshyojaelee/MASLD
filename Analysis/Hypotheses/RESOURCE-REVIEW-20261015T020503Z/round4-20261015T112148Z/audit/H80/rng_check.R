# Audit H80: which pipeline functions draw random numbers (susieR 0.14.2, coloc 5.2.3)
pat <- "sample[(]|runif|rnorm|set[.]seed"
for (fn in c("coloc.susie","coloc.bf_bf","coloc.signals")) {
  f <- deparse(get(fn, envir = asNamespace("coloc"))); h <- grep(pat, f, value = TRUE)
  cat("coloc::", fn, " random-draw lines: ", length(h), "\n", sep = ""); if (length(h)) cat(h, sep = "\n")
}
f <- deparse(susieR::susie_rss); cat("susie_rss -> ", grep("susie_suff_stat|susie_rss_lambda", f, value = TRUE), sep = "\n")
fs <- deparse(susieR::susie_suff_stat); cat("susie_suff_stat -> ", grep("susie_get_cs|n_purity", fs, value = TRUE), sep = "\n")
cat("get_purity draws sample(pos, n) when a candidate CS has more than n_purity =",
    formals(susieR::susie_get_cs)$n_purity, "variables\n")
