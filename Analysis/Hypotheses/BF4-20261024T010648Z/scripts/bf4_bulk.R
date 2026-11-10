#!/usr/bin/env Rscript
# BF4 step 2 (PRESPEC_BF4.md): cell-content marker scores in the 774 bulk biopsies; within-stage share of continuum
# position explained by bile-duct and fibroblast content; matched random-gene control; two-axis split; gene meta.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE); STEP1 <- args[1]; OUT <- args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
EXPR <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261009T171017Z",
                  "unsupervised/normalized_expression_by_cohort")
CG0 <- file.path(ROOT, "Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z")
SEED <- 20261024; NBOOT <- 2000; NNULL <- 1000; NMARK <- 30
COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
TYPES <- c(bile_duct = "Cholangiocytes", fibroblast = "Fibroblasts", endothelial = "Endothelial cells",
           macrophage = "Macrophages", t_cell = "T cells", hepatocyte = "Hepatocytes")
set.seed(SEED)
RES <- list()
istrue <- function(x) as.character(x) %in% c("TRUE", "True", "true")

u <- fread(file.path(CG0, "units.tsv"))[, .(sample_id, dataset, inferred_sex, fibrosis_stage, nas_score, axis_raw)]
genes <- fread(file.path(CG0, "genes.tsv"))[, .(g, gene_name, AveExpr, cont_z, signature_gene)]
genes[, signature_gene := istrue(signature_gene)]

# ---------------------------------------------------------------- expression, z scored within cohort
Z <- list(); RK <- list()
for (cc in COHORTS) {
  X <- readRDS(file.path(EXPR, paste0(cc, ".rds")))
  ids <- u[dataset == cc]$sample_id
  stopifnot(all(ids %in% colnames(X)))
  X <- X[, ids]
  X <- X[apply(X, 1, sd) > 0, ]
  Z[[cc]] <- t(scale(t(X)))
  RK[[cc]] <- t(apply(X, 1, rank))
  cat(cc, nrow(X), "genes", ncol(X), "biopsies\n")
}
common <- Reduce(intersect, lapply(Z, rownames))
cat("genes in all cohorts:", length(common), "\n")

# ---------------------------------------------------------------- marker sets
cand <- fread(file.path(STEP1, "marker_candidates.tsv"))
cand[, top10pct_abs_cont_z := istrue(top10pct_abs_cont_z)]
pick <- function(strict) {
  m <- cand[ensembl %in% common & (!strict | !top10pct_abs_cont_z)]
  m[order(rank), head(.SD, NMARK), by = cell_type]
}
MK <- list(main = pick(FALSE), s2 = pick(TRUE))
for (k in names(MK)) fwrite(MK[[k]], file.path(OUT, paste0("markers_", k, ".tsv")), sep = "\t")
RES$markers <- lapply(MK, function(m) m[, .(n = .N, genes = paste(gene, collapse = ",")), by = cell_type])

score <- function(cc, ens) { s <- colMeans(Z[[cc]][ens, , drop = FALSE]); as.numeric(scale(s)) }
D <- rbindlist(lapply(COHORTS, function(cc) {
  d <- u[dataset == cc]
  d[, y := as.numeric(scale(axis_raw))]
  for (k in names(TYPES)) {
    d[, (k) := score(cc, MK$main[cell_type == TYPES[[k]]]$ensembl)]
    d[, (paste0(k, "_s2")) := score(cc, MK$s2[cell_type == TYPES[[k]]]$ensembl)]
  }
  d
}))

# ---------------------------------------------------------------- P1: within-stage share explained
staged <- D[!is.na(fibrosis_stage)]
SC <- unique(staged$dataset)
partial <- function(d, vars) {
  per <- sapply(SC, function(cc) {
    x <- d[dataset == cc]
    a <- lm(y ~ factor(fibrosis_stage) + inferred_sex, data = x)
    b <- lm(as.formula(paste("y ~ factor(fibrosis_stage) + inferred_sex +", paste(vars, collapse = " + "))), data = x)
    c(rss_a = sum(resid(a)^2), rss_b = sum(resid(b)^2))
  })
  list(pooled = sum(per["rss_a", ] - per["rss_b", ]) / sum(per["rss_a", ]),
       per_cohort = as.list((per["rss_a", ] - per["rss_b", ]) / per["rss_a", ]))
}
boot <- function(vars) {
  replicate(NBOOT, {
    b <- staged[, .SD[sample.int(.N, .N, replace = TRUE)], by = dataset]
    partial(b, vars)$pooled
  })
}
run <- function(vars) {
  pt <- partial(staged, vars); bs <- boot(vars)
  list(pooled = pt$pooled, ci = unname(quantile(bs, c(0.025, 0.975))), per_cohort = pt$per_cohort,
       n = staged[, .N, by = dataset][, setNames(as.list(N), dataset)])
}
RES$P1 <- run(c("bile_duct", "fibroblast"))
RES$P1$reading <- with(RES$P1, if (pooled >= 0.5) "mostly portal-tract content" else if (pooled >= 0.2) "partly" else "little")
cat("P1", RES$P1$pooled, RES$P1$ci, "\n")

# ---------------------------------------------------------------- N1: expression-matched random gene sets
genes[, dec := cut(rank(AveExpr, ties.method = "first"), 10, labels = FALSE)]
pool <- genes[g %in% common & !signature_gene & !(g %in% unlist(lapply(MK, `[[`, "ensembl")))]
draw_like <- function(ens) {
  need <- genes[g %in% ens, .N, by = dec]
  unlist(lapply(seq_len(nrow(need)), function(i) sample(pool[dec == need$dec[i]]$g, need$N[i])))
}
null <- replicate(NNULL, {
  rb <- draw_like(MK$main[cell_type == "Cholangiocytes"]$ensembl)
  rf <- draw_like(MK$main[cell_type == "Fibroblasts"]$ensembl)
  d <- copy(D)  # score every biopsy in the cohort, as for the marker scores, then keep staged biopsies
  for (cc in SC) d[dataset == cc, `:=`(r1 = score(cc, rb), r2 = score(cc, rf))]
  partial(d[!is.na(fibrosis_stage)], c("r1", "r2"))$pooled
})
RES$N1 <- list(null_median = median(null), null_q95 = unname(quantile(null, 0.95)), null_max = max(null),
               observed_percentile = mean(null < RES$P1$pooled))
fwrite(data.table(null_partial_r2 = null), file.path(OUT, "n1_null.tsv"), sep = "\t")
cat("N1", unlist(RES$N1), "\n")

# ---------------------------------------------------------------- S1, S2
RES$S1_all_six <- run(names(TYPES))
RES$S2_strict_markers <- run(c("bile_duct_s2", "fibroblast_s2"))

# ---------------------------------------------------------------- S3: portal-content axis and residual axis
D[, c("P", "H") := {
  f <- lm(y ~ bile_duct + fibroblast); list(fitted(f), resid(f))
}, by = dataset]
fz <- function(r, n) list(r = r, n = n)
pool_r <- function(rs, ns) { w <- ns - 3; tanh(sum(w * atanh(rs)) / sum(w)) }
s3 <- list()
for (out in c("fibrosis_stage", "nas_score")) {
  tab <- D[!is.na(get(out)), .(n = .N, r_position = cor(y, get(out), method = "spearman"),
                               r_P = cor(P, get(out), method = "spearman"), r_H = cor(H, get(out), method = "spearman")),
           by = dataset]
  s3[[out]] <- list(per_cohort = tab, pooled = tab[, .(r_position = pool_r(r_position, n), r_P = pool_r(r_P, n),
                                                        r_H = pool_r(r_H, n), n = sum(n))])
}
s3$r2_position_on_P <- D[, .(r2 = summary(lm(y ~ bile_duct + fibroblast))$r.squared, n = .N), by = dataset]
RES$S3 <- s3
fwrite(D, file.path(OUT, "biopsy_scores.tsv"), sep = "\t")

# ---------------------------------------------------------------- S4: genes along each axis
meta <- function(axis) {
  zs <- sapply(COHORTS, function(cc) {
    v <- rank(D[dataset == cc][[axis]]); n <- length(v)
    r <- as.numeric(cor(t(RK[[cc]][common, ]), v))
    atanh(pmin(pmax(r, -0.999999), 0.999999)) * sqrt(n - 3)
  })
  ns <- D[, .N, by = dataset][match(COHORTS, dataset)]$N
  zm <- as.numeric(zs %*% sqrt(ns)) / sqrt(sum(ns))
  zm
}
gt <- data.table(g = common, z_H = meta("H"), z_P = meta("P"))
gt[, q_H := p.adjust(2 * pnorm(-abs(z_H)), "BH")][, q_P := p.adjust(2 * pnorm(-abs(z_P)), "BH")]
src <- fread(file.path(STEP1, "gene_top_celltype.tsv"))
gt <- merge(gt, genes[, .(g, gene_name, cont_z, signature_gene)], by = "g", all.x = TRUE)
gt <- merge(gt, src[, .(g = ensembl, top_type, hep_ratio)], by = "g", all.x = TRUE)
fwrite(gt[order(-z_H)], file.path(OUT, "gene_axes.tsv"), sep = "\t")
top <- function(col, dir, k) head(gt[order(dir * -get(col))], k)
RES$S4 <- list(
  n_genes = nrow(gt), n_q05_H = gt[q_H < 0.05, .N], n_q05_P = gt[q_P < 0.05, .N],
  H_up = top("z_H", 1, 25)$gene_name, H_down = top("z_H", -1, 25)$gene_name,
  P_up = top("z_P", 1, 25)$gene_name, P_down = top("z_P", -1, 25)$gene_name,
  named = gt[gene_name %in% c("AKR1B10", "CA12", "CIDEC", "ALDH3A1"), .(gene_name, z_H, q_H, z_P, q_P)],
  signature_genes = gt[signature_gene == TRUE, .(n = .N, median_z_H = median(z_H), median_z_P = median(z_P),
                                                 n_H_q05_pos = sum(q_H < 0.05 & z_H > 0), n_P_q05_pos = sum(q_P < 0.05 & z_P > 0))])

# ---------------------------------------------------------------- S5: which cells make the top genes of each axis
bg <- gt[!is.na(top_type), .N, by = top_type][, bg_share := N / sum(N)]
s5 <- lapply(c(H = "z_H", P = "z_P"), function(col) {
  t100 <- head(gt[!is.na(top_type)][order(-get(col))], 100)
  merge(t100[, .N, by = top_type][, share := N / 100], bg[, .(top_type, bg_share)], by = "top_type", all = TRUE)
})
RES$S5 <- s5

write_json(RES, file.path(OUT, "results.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(RES[c("P1", "N1", "S3", "S4", "S5")])
