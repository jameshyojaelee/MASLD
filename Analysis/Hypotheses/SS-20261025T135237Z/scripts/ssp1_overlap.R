# SS SSP1 (PRESPEC_SS_addendum_SSP1_overlap.md): DEG-list overlap vs sample size within SteatoSITE.
# Usage:
#   Rscript ssp1_overlap.R draws <out>              # make all draws (seed 20261027), write draws.tsv.gz
#   Rscript ssp1_overlap.R fit <out> <i> <k>        # fit subsets whose row index mod k == i (0-based), one rds each
#   Rscript ssp1_overlap.R summarize <out>          # per-draw statistics and summary
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R")
a <- commandArgs(trailingOnly = TRUE); mode <- a[1]; OUT <- a[2]
S2 <- file.path(SS, "out_s2_qc")
STAGE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv"
SIZES <- c(30L, 65L, 109L, 197L, 246L); B <- 50L

load_units <- function() {
  units <- fread(file.path(S2, "units.tsv")); st <- read_stage(STAGE)
  u <- merge(units[include == TRUE], st, by = "sample_alias", all.x = TRUE)
  u <- u[!is.na(stage_code) & !is.na(inferred_sex)]; setorder(u, sample_alias)
  stopifnot(nrow(u) == 493L); add_labels(u, u$stage_code)
}

if (mode == "draws") {
  OUT <- new_outdir(OUT); u <- load_units()
  set.seed(20261027L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  ok <- function(ids) { b <- u$bin[match(ids, u$sample_alias)]; sum(b == "F0F1") >= 3 && sum(b == "F3F4") >= 3 }
  rows <- list(); replaced <- 0L
  for (n in SIZES) for (d in seq_len(B)) repeat {
    o <- sample(u$sample_alias); A <- o[seq_len(n)]; Bs <- o[n + seq_len(n)]
    if (ok(A) && ok(Bs)) { rows[[length(rows) + 1L]] <- data.table(design = "pair", size = n, draw = d, subset = "A", sample_alias = A)
      rows[[length(rows) + 1L]] <- data.table(design = "pair", size = n, draw = d, subset = "B", sample_alias = Bs); break }
    replaced <- replaced + 1L }
  for (d in seq_len(B)) { o <- sample(u$sample_alias); g <- rep(1:5, length.out = length(o))
    for (k in 1:5) rows[[length(rows) + 1L]] <- data.table(design = "partition", size = sum(g == k), draw = d, subset = as.character(k), sample_alias = o[g == k]) }
  D <- rbindlist(rows); fwrite(D, file.path(OUT, "draws.tsv.gz"), sep = "\t")
  writeLines(c(paste("replaced pair draws:", replaced), paste("subsets:", uniqueN(D[, .(design, size, draw, subset)]))), file.path(OUT, "draw_notes.txt"))
  log_msg("draws written; replaced", replaced); quit(save = "no")
}

if (mode == "fit") {
  i <- as.integer(a[3]); k <- as.integer(a[4]); u <- load_units()
  dge <- readRDS(file.path(S2, "dge_included.rds"))
  if (!"group" %in% names(dge$samples)) dge$samples <- cbind(group = factor(rep(1L, ncol(dge))), dge$samples)
  D <- fread(file.path(OUT, "draws.tsv.gz")); keys <- unique(D[, .(design, size, draw, subset)])
  todo <- keys[(seq_len(nrow(keys)) - 1L) %% k == i]
  dir.create(file.path(OUT, "fits"), showWarnings = FALSE)
  res <- mclapply(seq_len(nrow(todo)), function(r) {
    kk <- todo[r]; f <- file.path(OUT, "fits", sprintf("%s_n%03d_d%02d_%s.rds", kk$design, kk$size, kk$draw, kk$subset))
    if (file.exists(f)) return(TRUE)
    ids <- D[design == kk$design & size == kk$size & draw == kk$draw & subset == kk$subset, sample_alias]
    info <- u[match(ids, sample_alias)]; d <- dge[, ids]
    keep <- filterByExpr(d, group = info$bin); d <- calcNormFactors(d[keep, , keep.lib.sizes = TRUE], method = "RLE")
    fit <- fit_qw(d, design_of(info, "bin"))
    tt <- topTable(fit, coef = "binF3F4", number = Inf, sort.by = "none")
    saveRDS(data.table(gene = rownames(tt), logFC = tt$logFC, q = tt$adj.P.Val), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
    TRUE }, mc.cores = NCORES, mc.preschedule = FALSE)
  stopifnot(all(vapply(res, isTRUE, TRUE))); log_msg("fit worker", i, "of", k, "done:", nrow(todo)); quit(save = "no")
}

if (mode == "summarize") {
  D <- fread(file.path(OUT, "draws.tsv.gz")); keys <- unique(D[, .(design, size, draw, subset)])
  rd <- function(kk) readRDS(file.path(OUT, "fits", sprintf("%s_n%03d_d%02d_%s.rds", kk$design, kk$size, kk$draw, kk$subset)))
  P <- rbindlist(lapply(SIZES, function(n) rbindlist(lapply(seq_len(B), function(d) {
    A <- rd(list(design = "pair", size = n, draw = d, subset = "A")); Bt <- rd(list(design = "pair", size = n, draw = d, subset = "B"))
    ca <- A[q < 0.05]; cb <- Bt[q < 0.05]; sh <- merge(ca, cb, by = "gene"); both <- merge(A, Bt, by = "gene")
    data.table(size = n, draw = d, calls_A = nrow(ca), calls_B = nrow(cb), shared = nrow(sh),
               jaccard = nrow(sh) / max(1, length(union(ca$gene, cb$gene))),
               frac_called_in_other = mean(c(nrow(sh) / max(1, nrow(ca)), nrow(sh) / max(1, nrow(cb)))),
               sign_discordant_shared = if (nrow(sh)) mean(sign(sh$logFC.x) != sign(sh$logFC.y)) else NA_real_,
               spearman_logfc = cor(both$logFC.x, both$logFC.y, method = "spearman")) }))))
  Q <- rbindlist(lapply(seq_len(B), function(d) {
    C <- rbindlist(lapply(1:5, function(k) rd(list(design = "partition", size = keys[design == "partition" & draw == d & subset == k, size], draw = d, subset = as.character(k)))[q < 0.05][, sub := k]))
    g <- C[, .(n_sub = .N, signs = uniqueN(sign(logFC))), by = gene]
    data.table(draw = d, union = nrow(g), frac_single = mean(g$n_sub == 1L), multi = sum(g$n_sub >= 2L),
               sign_discordant_multi = if (any(g$n_sub >= 2L)) mean(g[n_sub >= 2L]$signs > 1L) else NA_real_) }))
  s <- function(x) c(median = median(x, na.rm = TRUE), lo = unname(quantile(x, 0.025, na.rm = TRUE)), hi = unname(quantile(x, 0.975, na.rm = TRUE)))
  SP <- P[, as.list(unlist(lapply(.SD, s))), by = size, .SDcols = setdiff(names(P), c("size", "draw"))]
  SQ <- as.data.table(as.list(unlist(lapply(Q[, -"draw"], s))))
  fwrite(P, file.path(OUT, "pairs_per_draw.tsv"), sep = "\t"); fwrite(Q, file.path(OUT, "partitions_per_draw.tsv"), sep = "\t")
  fwrite(SP, file.path(OUT, "pairs_summary.tsv"), sep = "\t"); fwrite(SQ, file.path(OUT, "partitions_summary.tsv"), sep = "\t")
  writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  print(SP); print(SQ); log_msg("summary written")
}
