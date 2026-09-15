#!/usr/bin/env Rscript
# Null B (sensitivity): caliper / nearest-neighbour matching on standardized
# discovery mean and log-variance. Emits gene identifiers only. Fits nothing.
suppressPackageStartupMessages({library(digest)})
UNIV <- Sys.getenv("MNC_UNIVERSE_DIR"); OUT <- Sys.getenv("MNC_OUT_DIR")
stopifnot(nzchar(UNIV), nzchar(OUT))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
stopifnot(!file.exists(file.path(OUT, "null_gene_sets_caliper.tsv.gz")))
SEED_BASE <- 20270828L; N_DRAWS <- 1000L; K_BASE <- 50L
say <- function(...) cat(sprintf(...), "\n", sep = "")

U  <- read.delim(file.path(UNIV, "sampling_universe.tsv"), stringsAsFactors = FALSE)
S  <- read.delim(file.path(UNIV, "signature_strata.tsv"),  stringsAsFactors = FALSE)
stopifnot(nrow(S) == 139L, nrow(U) == 13810L, !any(S$gene_id_base %in% U$gene_id_base))

mu_m <- mean(U$discovery_mean); mu_s <- stats::sd(U$discovery_mean)
lv_u <- log(U$discovery_var);  lv_m <- mean(lv_u); lv_s <- stats::sd(lv_u)
zU <- cbind((U$discovery_mean - mu_m) / mu_s, (lv_u - lv_m) / lv_s)
zS <- cbind((S$discovery_mean - mu_m) / mu_s, (log(S$discovery_var) - lv_m) / lv_s)

ord <- order(S$gene_id_base); S <- S[ord, ]; zS <- zS[ord, , drop = FALSE]
nn <- lapply(seq_len(nrow(zS)), function(i) {
  d <- sqrt((zU[, 1] - zS[i, 1])^2 + (zU[, 2] - zS[i, 2])^2)
  order(d)[seq_len(200L)]
})
cal <- vapply(seq_len(nrow(zS)), function(i) {
  d <- sqrt((zU[nn[[i]][seq_len(K_BASE)], 1] - zS[i, 1])^2 +
            (zU[nn[[i]][seq_len(K_BASE)], 2] - zS[i, 2])^2); max(d)
}, numeric(1))
say("caliper radius over 139 anchors: median %.4f  max %.4f", median(cal), max(cal))

draw_one <- function(b) {
  set.seed(SEED_BASE + b); picked <- integer(0)
  for (i in seq_len(139L)) {
    got <- NA_integer_
    for (K in c(K_BASE, 100L, 200L)) {
      cand <- setdiff(nn[[i]][seq_len(K)], picked)
      if (length(cand)) { got <- if (length(cand) == 1L) cand else sample(cand, 1L); break }
    }
    if (is.na(got)) stop(sprintf("draw %d anchor %d exhausted 200 neighbours", b, i))
    picked <- c(picked, got)
  }
  stopifnot(length(picked) == 139L, !anyDuplicated(picked))
  U$gene_id_base[picked]
}
sets <- lapply(seq_len(N_DRAWS), draw_one)
per_sha <- vapply(sets, function(s) digest(paste(sort(s), collapse = "\n"),
                                           algo = "sha256", serialize = FALSE), character(1))
stopifnot(!anyDuplicated(per_sha))
gz <- gzfile(file.path(OUT, "null_gene_sets_caliper.tsv.gz"), "w")
write.table(data.frame(draw_index = rep(seq_len(N_DRAWS), each = 139L),
                       gene_id_base = unlist(sets), stringsAsFactors = FALSE),
            gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
write.table(data.frame(draw_index = seq_len(N_DRAWS), set_sha256 = per_sha,
                       stringsAsFactors = FALSE),
            file.path(OUT, "null_gene_set_digests_caliper.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(anchor_gene_id_base = S$gene_id_base, caliper_radius_k50 = cal,
                       stringsAsFactors = FALSE),
            file.path(OUT, "caliper_radii.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(
  key = c("seed_base", "n_draws", "k_base", "caliper_radius_median", "caliper_radius_max",
          "gene_sets_collection_sha256"),
  value = c(SEED_BASE, N_DRAWS, K_BASE, sprintf("%.6f", median(cal)), sprintf("%.6f", max(cal)),
            digest(paste(per_sha, collapse = "\n"), algo = "sha256", serialize = FALSE)),
  stringsAsFactors = FALSE),
  file.path(OUT, "caliper_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_caliper.txt"))
say("DONE")
