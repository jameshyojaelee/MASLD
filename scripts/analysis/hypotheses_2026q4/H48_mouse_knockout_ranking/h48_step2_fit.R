#!/usr/bin/env Rscript
# H48 step 2 (Amendment 2 fallback): per-knockout limma-voom effects on (sgRNA, letter)
# pseudobulks, design ~ letter + group, KO guides vs all NC guides; NC pseudo-knockout draws.
# Writes OUT/step2/. No scoring here.
suppressPackageStartupMessages({
  library(Matrix); library(edgeR); library(limma); library(data.table); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]; NCORES <- as.integer(args[2])
TEST <- length(args) >= 3 && args[3] == "test"
B <- if (TEST) 3L else 200L; B_FULL <- if (TEST) 2L else 20L; SEED <- 48L; MINCELLS <- 10L; DROP_LETTER <- "A"
s1 <- file.path(OUT, "step1"); s2 <- file.path(OUT, "step2")
dir.create(s2, showWarnings = FALSE)
if (length(list.files(s2))) stop("step2 directory not empty; refusing to overwrite")

univ <- fread(file.path(s1, "gene_universe.tsv"))
PB <- as(readMM(file.path(s1, "pseudobulk_counts.mtx")), "CsparseMatrix")
cols <- fread(file.path(s1, "pseudobulk_columns.tsv"))
stopifnot(nrow(PB) == nrow(univ), ncol(PB) == nrow(cols))
rownames(PB) <- univ$gene_id; colnames(PB) <- cols$key
NCM <- as(readMM(file.path(s1, "nc_cells_counts.mtx")), "CsparseMatrix")
ncm <- fread(file.path(s1, "nc_cells_meta.tsv"))
stopifnot(ncol(NCM) == nrow(ncm))
elig <- fread(file.path(s1, "knockout_eligibility.tsv"))
# tested set plus descriptive-only set (Amendment 3): >= 50 cells and exactly one guide with >= 10 cells
elig[, analysis := fifelse(tested, "tested", fifelse(n_cells_total >= 50 & n_guides_ge10 == 1, "descriptive_single_guide", "not_fit"))]
fitset <- elig[analysis != "not_fit"]
if (TEST) fitset <- rbind(head(fitset[analysis == "tested"], 3), head(fitset[analysis != "tested"], 1))
sigs <- rbindlist(lapply(list.files(s1, "^signature_.*tsv$", full.names = TRUE), fread), fill = TRUE)
sig_genes <- sort(unique(sigs$mouse_gene_id))
sig_idx <- match(sig_genes, univ$gene_id); stopifnot(!anyNA(sig_idx))
writeLines(sig_genes, file.path(s2, "sig_union_genes.txt"))

nc_cols <- cols[target_gene == "control" & included == TRUE]
NC <- as.matrix(PB[, nc_cols$key, drop = FALSE])
nc_letter <- nc_cols$letter
cells_by_letter <- split(seq_len(nrow(ncm)), ncm$letter)
nc_key_idx <- match(ncm$key, nc_cols$key); stopifnot(!anyNA(nc_key_idx))

fit_lfc <- function(ko_mat, ko_letter, ctl_mat, ctl_letter) {
  y <- cbind(ko_mat, ctl_mat)
  grp <- factor(c(rep("KO", ncol(ko_mat)), rep("NC", ncol(ctl_mat))), levels = c("NC", "KO"))
  le <- factor(c(ko_letter, ctl_letter))
  design <- if (nlevels(le) > 1) model.matrix(~ le + grp) else model.matrix(~ grp)
  d <- calcNormFactors(DGEList(y))
  v <- voom(d, design, plot = FALSE)
  fit <- lmFit(v, design)
  unname(fit$coefficients[, "grpKO"])
}

one_ko <- function(i) {
  g <- fitset$target_gene[i]
  kc <- cols[target_gene == g & included == TRUE]
  if (fitset$analysis[i] == "descriptive_single_guide") {
    gc <- kc[, .(n = sum(n_cells)), by = sgrna_id]
    kc <- kc[sgrna_id %in% gc[n >= MINCELLS, sgrna_id]]
  }
  ko <- as.matrix(PB[, kc$key, drop = FALSE])
  real <- fit_lfc(ko, kc$letter, NC, nc_letter)
  ka <- kc$letter != DROP_LETTER; na <- nc_letter != DROP_LETTER
  real_dropA <- if (sum(ka) >= 2) fit_lfc(ko[, ka, drop = FALSE], kc$letter[ka], NC[, na, drop = FALSE], nc_letter[na]) else rep(NA_real_, nrow(PB))
  null_full <- matrix(NA_real_, nrow(PB), B_FULL)
  null_sig <- matrix(NA_real_, length(sig_idx), B)
  if (fitset$analysis[i] == "tested") {
    for (b in seq_len(B)) {
      set.seed(SEED + 1000L * i + b)
      pk <- matrix(0, nrow(PB), nrow(kc)); ctl <- NC
      for (L in unique(kc$letter)) {
        js <- which(kc$letter == L)
        pool <- cells_by_letter[[L]]
        take <- pool[sample.int(length(pool), sum(kc$n_cells[js]))]
        split_at <- rep(js, kc$n_cells[js])
        for (j in js) {
          idx <- take[split_at == j]
          pk[, j] <- Matrix::rowSums(NCM[, idx, drop = FALSE])
        }
        sub <- as.matrix(NCM[, take, drop = FALSE] %*% sparseMatrix(i = seq_along(take), j = nc_key_idx[take], x = 1,
                                                                      dims = c(length(take), ncol(NC))))
        ctl <- ctl - sub
      }
      keepc <- colSums(ctl) > 0
      l <- fit_lfc(pk, kc$letter, ctl[, keepc, drop = FALSE], nc_letter[keepc])
      null_sig[, b] <- l[sig_idx]
      if (b <= B_FULL) null_full[, b] <- l
    }
  }
  list(gene = g, analysis = fitset$analysis[i], n_samples = nrow(kc), real = real, real_dropA = real_dropA,
       null_full = null_full, null_sig = null_sig)
}

t0 <- Sys.time()
res <- mclapply(seq_len(nrow(fitset)), function(i) tryCatch(one_ko(i), error = function(e) list(err = conditionMessage(e), i = i)),
                mc.cores = NCORES, mc.preschedule = FALSE)
bad <- vapply(res, function(r) !is.null(r$err) || is.null(r$real), logical(1))
if (any(bad)) { print(res[bad]); stop("failed fits: ", sum(bad)) }
message("fits done in ", format(Sys.time() - t0))
genes <- vapply(res, `[[`, "", "gene")
eff <- do.call(cbind, lapply(res, `[[`, "real")); colnames(eff) <- genes
effa <- do.call(cbind, lapply(res, `[[`, "real_dropA")); colnames(effa) <- genes
fwrite(data.table(gene_id = univ$gene_id, eff), file.path(s2, "effects_real_log2fc.tsv.gz"), sep = "\t")
fwrite(data.table(gene_id = univ$gene_id, effa), file.path(s2, "effects_real_dropA_log2fc.tsv.gz"), sep = "\t")
fwrite(data.table(target_gene = genes, analysis = vapply(res, `[[`, "", "analysis"),
                  n_ko_samples = vapply(res, `[[`, 0L, "n_samples")), file.path(s2, "fit_set.tsv"), sep = "\t")
wr <- function(lst, field, path) {  # per KO: genes x draws column-major = C order [draw, gene]
  con <- file(path, "wb")
  for (r in lst) writeBin(as.numeric(r[[field]]), con, size = 4)
  close(con)
}
wr(res, "null_full", file.path(s2, "null_full_f32.bin"))
wr(res, "null_sig", file.path(s2, "null_sig_f32.bin"))
jsonlite::write_json(list(kos = genes, n_genes = nrow(PB), n_sig_union = length(sig_idx), B = B, B_full = B_FULL,
                          drop_letter = DROP_LETTER, seed = SEED, n_nc_samples = ncol(NC)),
                     file.path(s2, "step2_dims.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(s2, "sessionInfo.txt"))
message("done")
