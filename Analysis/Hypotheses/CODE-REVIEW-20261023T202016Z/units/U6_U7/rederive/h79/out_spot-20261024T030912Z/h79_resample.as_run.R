#!/usr/bin/env Rscript
# Code review U7 (H79): independent re-derivation of the T1 and T2 reference distributions.
# Written from METHODS.md "Genetic ancestry sensitivity", not from the producer script.
#   T1: ~ dataset + inferred_sex + PC1 + PC2 + PC3 + PC4 + group_binary; reference = PCs permuted jointly
#       within the 10 cohort-by-group cells (B1 = 200).
#   T2: ~ dataset + inferred_sex + group_binary on the 614 EUR units; reference = random within-cell subsets
#       of the EUR cell sizes (B2 = 100).
#   Each fit: voomWithQualityWeights, lmFit, eBayes; call = BH adj p < 0.05 and |log2FC| > 0.5 over the
#   24,115-gene universe; statistic C = number of genes whose call differs from the stored 1,854 calls.
#   p = (1 + k)/(B + 1), k = reference draws with C >= observed; Holm over {T1, T2}.
# Modes:
#   spot  seed 20260923: compare draw indices with the producer's draw_indices.rds; refit the true T1 and T2
#         fits, reference draws 1-4 (T1) and 1-2 (T2), and T2 with RLE factors recomputed on the EUR subset.
#   full  any seed: all B1 + B2 reference fits, k, p, Holm, the ancestry-sensitive rule (<= 9 of 200).
# Observed C values (233 and 570) are taken from the spot run, not from the producer.
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
ga <- function(k, d = NULL) { i <- match(paste0("--", k), args); if (is.na(i)) d else args[i + 1L] }
SEED <- as.integer(ga("seed")); OUT <- ga("out"); CORES <- as.integer(ga("cores", "1"))
MODE <- ga("mode", "full"); B1 <- as.integer(ga("b1", "200")); B2 <- as.integer(ga("b2", "100"))
OBS_C1 <- as.integer(ga("obs-c1", NA)); OBS_C2 <- as.integer(ga("obs-c2", NA))
stopifnot(!is.na(SEED), !is.null(OUT))
dir.create(file.path(OUT, "draws"), recursive = TRUE, showWarnings = FALSE)
say <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")), ..., "\n", sep = ""); flush(stdout()) }

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
FIT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANC <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv")
PROD_DRAWS <- file.path(ROOT, "Analysis/Hypotheses/H79-20261015T125125Z/results/run-20261015T132031Z/draw_indices.rds")
in_files <- c(dge = file.path(FIT, "merged_dge.rds"), meta = file.path(FIT, "meta_matched.rds"),
              deg = file.path(FIT, "deg_results.csv"), anc = ANC)
fwrite(data.table(role = names(in_files), path = in_files,
                  sha256 = vapply(in_files, function(p) sub(" .*", "", system2("sha256sum", shQuote(p), stdout = TRUE)), "")),
       file.path(OUT, "inputs_sha256.tsv"), sep = "\t")

dge <- readRDS(in_files[["dge"]])
meta <- as.data.table(readRDS(in_files[["meta"]]))
deg <- fread(in_files[["deg"]])
anc <- fread(ANC)
N <- ncol(dge); G <- nrow(dge)
say("DGEList ", G, " genes x ", N, " units")

# Join by identifier (log counts; no inner-join drops allowed)
im <- match(colnames(dge), meta$sample_id); ia <- match(colnames(dge), anc$run); ig <- match(rownames(dge), deg$gene)
say("join: meta matched ", sum(!is.na(im)), "/", N, "; ancestry matched ", sum(!is.na(ia)), "/", N,
    "; deg rows matched ", sum(!is.na(ig)), "/", G)
stopifnot(!anyNA(im), !anyNA(ia), !anyNA(ig))
meta <- meta[im]; anc <- anc[ia]; deg <- deg[ig]
stopifnot(all(meta$group_binary == dge$samples$group_binary), all(meta$dataset == dge$samples$dataset))
PCN <- paste0("PC", 1:4)
d0 <- data.frame(dataset = factor(meta$dataset), inferred_sex = factor(meta$inferred_sex),
                 group_binary = factor(meta$group_binary, levels = c("Control", "Disease")))
stopifnot(!anyNA(d0), !anyNA(as.matrix(anc[, ..PCN])))
PC <- as.matrix(anc[, ..PCN])
cell <- interaction(d0$dataset, d0$group_binary, lex.order = TRUE, drop = TRUE)
eur <- anc$ancestry_call == "EUR"
say("units: control ", sum(d0$group_binary == "Control"), ", disease ", sum(d0$group_binary == "Disease"),
    "; EUR ", sum(eur), " (control ", sum(eur & d0$group_binary == "Control"), ", disease ", sum(eur & d0$group_binary == "Disease"), ")")
adopted <- !is.na(deg$padj) & deg$padj < 0.05 & abs(deg$logFC) > 0.5
say("stored calls ", sum(adopted), " (up ", sum(adopted & deg$logFC > 0), ", down ", sum(adopted & deg$logFC < 0), ")")

# Draws: one seed, all T1 permutations first, then all T2 subsets, before any fit
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
cells <- split(seq_len(N), cell)
perms <- vector("list", B1)
for (b in seq_len(B1)) {
  p <- seq_len(N)
  for (ix in cells) p[ix] <- ix[sample.int(length(ix))]
  perms[[b]] <- p
}
n_eur_cell <- vapply(cells, function(ix) sum(eur[ix]), 0L)
subs <- vector("list", B2)
for (b in seq_len(B2)) {
  s <- integer(0)
  for (cn in names(cells)) { ix <- cells[[cn]]; s <- c(s, ix[sample.int(length(ix), n_eur_cell[[cn]])]) }
  subs[[b]] <- sort(s)
}
saveRDS(list(seed = SEED, perms = perms, subs = subs, cells = names(cells)), file.path(OUT, "draws_mine.rds"))
if (MODE %in% c("spot", "draws")) {
  pd <- readRDS(PROD_DRAWS)
  cmp <- data.table(check = c("perm_identical_all", "subs_identical_all", "n_perm_prod", "n_subs_prod"),
                    value = c(identical(lapply(pd$perm, as.integer), lapply(perms, as.integer)),
                              identical(lapply(pd$subs, as.integer), lapply(subs, as.integer)),
                              length(pd$perm), length(pd$subs)))
  fwrite(cmp, file.path(OUT, "draws_vs_producer.tsv"), sep = "\t")
  say("draws vs producer: ", paste(cmp$check, cmp$value, sep = "=", collapse = "; "))
  if (MODE == "draws") quit(save = "no")
}

fit_calls <- function(dg, d, form) {
  X <- model.matrix(form, data = d)
  stopifnot(qr(X)$rank == ncol(X))
  v <- voomWithQualityWeights(dg, X)
  fit <- eBayes(lmFit(v, X))
  tt <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  stopifnot(identical(rownames(tt), rownames(dg)))
  list(logFC = tt$logFC, padj = tt$adj.P.Val)
}
F_T1 <- ~ dataset + inferred_sex + PC1 + PC2 + PC3 + PC4 + group_binary
F_T2 <- ~ dataset + inferred_sex + group_binary
summ <- function(id, f, n_units) {
  cl <- !is.na(f$padj) & f$padj < 0.05 & abs(f$logFC) > 0.5
  list(id = id, n_units = n_units, calls = sum(cl), lost = sum(adopted & !cl), gained = sum(!adopted & cl),
       C = sum(adopted != cl), changed_idx = which(adopted != cl), n_na_padj = sum(is.na(f$padj)))
}
run_task <- function(tk) {
  fn <- file.path(OUT, "draws", paste0(tk$id, ".rds"))
  if (file.exists(fn)) return("exists")
  t0 <- Sys.time()
  if (tk$kind == "t1") {
    d <- d0; d[PCN] <- as.data.frame(PC[tk$idx, , drop = FALSE]); r <- summ(tk$id, fit_calls(dge, d, F_T1), N)
  } else if (tk$kind == "t2") {
    u <- tk$idx; r <- summ(tk$id, fit_calls(dge[, u, keep.lib.sizes = TRUE], d0[u, ], F_T2), length(u))
  } else if (tk$kind == "t2_rle") {
    u <- tk$idx; dg <- dge[, u, keep.lib.sizes = TRUE]
    dg <- calcNormFactors(dg, method = "RLE")
    r <- summ(tk$id, fit_calls(dg, d0[u, ], F_T2), length(u))
  }
  r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(r, paste0(fn, ".part")); file.rename(paste0(fn, ".part"), fn)
  "done"
}
tasks <- list()
if (MODE == "spot") {
  tasks <- c(list(list(id = "t1_true", kind = "t1", idx = seq_len(N)),
                  list(id = "t2_true", kind = "t2", idx = which(eur)),
                  list(id = "t2_true_rle_recomputed", kind = "t2_rle", idx = which(eur))),
             lapply(1:4, function(b) list(id = sprintf("perm_%03d", b), kind = "t1", idx = perms[[b]])),
             lapply(1:2, function(b) list(id = sprintf("sub_%03d", b), kind = "t2", idx = subs[[b]])))
} else {
  tasks <- c(lapply(seq_len(B1), function(b) list(id = sprintf("perm_%03d", b), kind = "t1", idx = perms[[b]])),
             lapply(seq_len(B2), function(b) list(id = sprintf("sub_%03d", b), kind = "t2", idx = subs[[b]])))
}
say(MODE, ": ", length(tasks), " fits on ", CORES, " cores; seed ", SEED, " B1 ", B1, " B2 ", B2)
res <- mclapply(tasks, function(tk) tryCatch(run_task(tk), error = function(e) paste("ERROR", tk$id, conditionMessage(e))),
                mc.cores = CORES, mc.preschedule = FALSE)
bad <- grep("^ERROR", unlist(res), value = TRUE)
if (length(bad)) { say("errors: ", paste(head(bad, 5), collapse = " | ")); quit(status = 1) }
rd <- function(id) readRDS(file.path(OUT, "draws", paste0(id, ".rds")))
tab <- rbindlist(lapply(tasks, function(tk) { r <- rd(tk$id); as.data.table(r[c("id", "n_units", "calls", "lost", "gained", "C", "n_na_padj", "seconds")]) }))
fwrite(tab, file.path(OUT, "fit_summary.tsv"), sep = "\t")

if (MODE == "full") {
  c1 <- tab[startsWith(id, "perm_"), C]; c2 <- tab[startsWith(id, "sub_"), C]
  k1 <- sum(c1 >= OBS_C1); k2 <- sum(c2 >= OBS_C2)
  p1 <- (1 + k1) / (B1 + 1); p2 <- (1 + k2) / (B2 + 1)
  ph <- p.adjust(c(p1, p2), "holm")
  # ancestry-sensitive rule: changed with true PCs and changed in <= 9 of the B1 permuted fits
  chg <- matrix(FALSE, G, B1)
  for (b in seq_len(B1)) chg[rd(sprintf("perm_%03d", b))$changed_idx, b] <- TRUE
  nperm <- rowSums(chg)
  true_t1 <- ga("true-t1-rds")
  d1 <- NA_integer_
  if (!is.null(true_t1) && file.exists(true_t1)) {
    tchg <- logical(G); tchg[readRDS(true_t1)$changed_idx] <- TRUE
    d1 <- sum(tchg & nperm <= 9)
  }
  ff <- vapply(seq_len(B1), function(b) sum(chg[, b] & (nperm - chg[, b]) <= 9), 0)
  out <- data.table(seed = SEED, B1 = B1, B2 = B2, obs_C1 = OBS_C1, obs_C2 = OBS_C2,
                    ref1_median = median(c1), ref1_q95 = unname(quantile(c1, 0.95)), ref1_max = max(c1), k1 = k1, p1 = p1,
                    ref2_median = median(c2), ref2_q95 = unname(quantile(c2, 0.95)), ref2_max = max(c2), k2 = k2, p2 = p2,
                    holm1 = ph[1], holm2 = ph[2], reject1 = ph[1] < 0.05, reject2 = ph[2] < 0.05,
                    d1_ancestry_sensitive = d1, false_flags_median = median(ff), false_flags_q95 = unname(quantile(ff, 0.95)))
  fwrite(out, file.path(OUT, "holm_result.tsv"), sep = "\t")
  say("RESULT ", paste(names(out), unlist(out), sep = "=", collapse = " "))
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("done")
