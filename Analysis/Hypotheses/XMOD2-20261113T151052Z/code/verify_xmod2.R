# XMOD2 checks: (1) trio regression vs CONF1, (2) independent recomputation of PPM1K and CFHR5 in GSE256398
# with lm() and an lm-based HC3, (3) GSE289173 betas vs CONF1 all_gene_betas, (4) donor counts vs trio.
# Usage: Rscript verify_xmod2.R <xmod2_out_dir> <check_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR)})
set.seed(20261109)
args <- commandArgs(trailingOnly = TRUE)
X <- args[1]; OUT <- args[2]
if (dir.exists(OUT)) stop("check dir exists"); dir.create(OUT)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CONF1 <- file.path(ROOT, "Analysis/Hypotheses/CONF1-20261109T221846Z")
V3 <- file.path(ROOT, "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk")
wr <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t")

maxdiff <- function(a, b, keys) {
  m <- merge(a, b, by = keys, suffixes = c(".new", ".old"))
  num <- intersect(names(a)[vapply(a, is.numeric, TRUE)], names(b)[vapply(b, is.numeric, TRUE)])
  num <- setdiff(num, keys)
  data.table(column = num, max_abs_diff = vapply(num, function(k) max(abs(m[[paste0(k, ".new")]] - m[[paste0(k, ".old")]]), na.rm = TRUE), 1),
             n_rows_new = nrow(a), n_rows_old = nrow(b), n_matched = nrow(m))
}
r1 <- maxdiff(fread(file.path(X, "trio_c2/per_cohort.tsv")), fread(file.path(CONF1, "out_c2_methionine/per_cohort.tsv")), c("cohort", "gene"))[, table := "c2 per_cohort"]
r2 <- maxdiff(fread(file.path(X, "trio_c5/primary.tsv")), fread(file.path(CONF1, "out_c5_trio/primary.tsv")), "gene")[, table := "c5 primary"]
r3 <- maxdiff(fread(file.path(X, "trio_c2/combined.tsv"))[, .(gene, analysis, T, p_perm, mde)],
              fread(file.path(CONF1, "out_c2_methionine/combined.tsv"))[, .(gene, analysis, T, p_perm, mde)], c("gene", "analysis"))[, table := "c2 combined"]
reg <- rbind(r1, r2, r3); print(reg); wr(reg, "regression_trio.tsv")

# independent recomputation, GSE256398
m <- fread(file.path(V3, "donor_metadata.tsv")); nn <- fread(file.path(V3, "nuclei_per_donor_lineage.tsv"))
pb <- fread(file.path(V3, "pseudobulk_hepatocyte.tsv.gz"))
m[, n_hep := nn$hepatocyte[match(donor, nn$donor)]]
el <- m[n_hep >= 30]
cn <- as.matrix(pb[, el$donor, with = FALSE]); rownames(cn) <- pb$symbol
y <- DGEList(cn); y <- calcNormFactors(y, method = "TMM")
lib <- y$samples$lib.size * y$samples$norm.factors
hc3_lm <- function(fit) { Xm <- model.matrix(fit); e <- residuals(fit); h <- hatvalues(fit)
  XtXi <- solve(crossprod(Xm)); V <- XtXi %*% crossprod(Xm * (e / (1 - h))) %*% XtXi; sqrt(diag(V)) }
ind <- rbindlist(lapply(c("PPM1K", "CFHR5"), function(g) {
  i <- which(pb$symbol == g); stopifnot(length(i) == 1)
  # log2 CPM with prior count 1 as in edgeR::cpm(log=TRUE): prior scaled by lib/mean(lib), lib augmented by 2*prior
  pc <- 1 * lib / mean(lib); lcpm <- log2((cn[i, ] + pc) / (lib + 2 * pc) * 1e6)
  d <- data.table(donor = el$donor, y = lcpm, stage = el$stage_ordinal, female = as.integer(el$sex == "female"), ln = log10(el$n_hep))
  d <- d[stage == 0 | stage >= 2][, dis := as.integer(stage >= 2)]
  f <- lm(y ~ dis + female + ln, d); se <- hc3_lm(f)
  data.table(cohort = "GSE256398", gene = g, n = nrow(d), n0 = sum(d$dis == 0), n1 = sum(d$dis == 1), beta_indep = coef(f)[["dis"]], se_indep = se[["dis"]])
}))
pc <- fread(file.path(X, "examples_c2/per_cohort.tsv"))
ind <- merge(ind, pc[, .(cohort, gene, n_script = n, beta_script = beta, se_script = se)], by = c("cohort", "gene"))
ind[, `:=`(beta_diff = beta_indep - beta_script, se_diff = se_indep - se_script)]
print(ind); wr(ind, "independent_recompute.tsv")

# GSE289173 vs CONF1 genome-wide betas
ab <- fread(file.path(CONF1, "out_c5_trio/all_gene_betas.tsv.gz"))
c5 <- fread(file.path(X, "examples_c5/primary.tsv"))
cmp <- merge(c5[, .(gene, beta_xmod2 = beta)], ab[, .(gene = symbol, beta_allgenes = beta)], by = "gene", all.x = TRUE)
cmp[, diff := beta_xmod2 - beta_allgenes]; print(cmp); wr(cmp, "gse289173_vs_all_gene_betas.tsv")

# donor counts per cohort: examples vs trio
dn <- merge(unique(fread(file.path(X, "examples_c2/per_cohort.tsv"))[, .(cohort, n, n0, n1)]),
            unique(fread(file.path(CONF1, "out_c2_methionine/per_cohort.tsv"))[, .(cohort, n, n0, n1)]), by = "cohort", suffixes = c(".ex", ".trio"), all = TRUE)
dn <- rbind(dn, data.table(cohort = "GSE289173", unique(c5[, .(n.ex = n, n0.ex = n0, n1.ex = n1)]),
                           unique(fread(file.path(CONF1, "out_c5_trio/primary.tsv"))[, .(n.trio = n, n0.trio = n0, n1.trio = n1)])), use.names = TRUE)
print(dn); wr(dn, "donor_counts.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
