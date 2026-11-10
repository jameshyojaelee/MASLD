# DIG1 job 67b (post hoc; rerun of the job 67 permutation part after its p_perm bug: the |obs| vector was recycled down
# columns, so every p was NA. z-scores of job 67 were read before this script). Job 67 header follows.
# DIG1 job 67 (post hoc, after job 66 was read). Is the genotype -> single-cell-state overlap specific, or a global
# shift in RNA composition?
# Read before this script: job 66 cameraPR. CCL20_STATE_UP (254 genes) rises with PNPLA3 148M (p 5e-10) and falls with
# the HSD17B13 junction allele (p 1e-11). But the HSD17B13 allele lowers 39 of 50 Hallmarks at BH 0.05 (OXPHOS, MYC,
# E2F, mTORC1 ...), and its top up genes are non-coding. That pattern can come from a library-composition axis.
# 1. Technical metrics per sample from STAR (Log.final.out, ReadsPerGene.out.tab): unique %, multi %, splices per
#    unique read, mismatch rate, fraction of unique reads assigned to genes (strand column with most assigned),
#    protein-coding share of release counts (GENCODE v49 gene_type). Metric ~ allele + dataset + sex (HC3);
#    INFLAMED_HEP ~ allele + dataset + sex + metrics.
# 2. Within-cohort permutation null for the gene-set statistics of the pooled genotype scan (log2 TMM CPM, prior 2,
#    genes by filterByExpr(group = dataset); ordinary t from expression ~ allele + dataset + sex via residualization).
#    Statistics: mean t of the set minus mean t of (a) all other genes, (b) expression-matched protein-coding genes
#    (20 AveExpr bins). Sets: CCL20_STATE_UP, CCL20_STATE_DN (job 66 definition), all Hallmark genes, protein-coding.
#    10,000 permutations of the allele within cohort, seed 20260923. Two-sided empirical p, (1 + #|null| >= |obs|)/(1 + B).
#    The same per cohort (allele + sex), 2,000 permutations each.
# 67b addition: is the genotype footprint on the state genes larger than its footprint on disease genes in general?
#    Gene-level genotype t is regressed on a fixed disease t vector (release deg_results t, disease vs control) and,
#    separately, on a fibrosis-stage t vector (expression ~ fibrosis + dataset + sex, same pipeline); statistic = mean
#    residual of the set minus expression-matched protein-coding genes; same within-cohort permutation null.
# Unit: one bulk sample per participant. No ancestry covariates. Cross-sectional.
# Usage: Rscript 67b_genotype_scan_null.R <job63_out_dir> <job66_out_dir> <gene_type_tsv> <new_out_dir>
suppressPackageStartupMessages({library(edgeR); library(data.table)})
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 4)
g63 <- args[1]; g66 <- args[2]; gtype <- args[3]; out <- args[4]; stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
set.seed(20260923)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
rel <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
GMT <- file.path(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/msigdb.v2025.1.Hs.symbols.gmt")
D <- read.delim(file.path(g63, "sample_scores_genotypes.tsv"))
d <- readRDS(file.path(rel, "merged_dge.rds")); stopifnot(identical(colnames(d$counts), as.character(D$sample_id)))
deg <- read.csv(file.path(rel, "deg_results.csv")); sym <- setNames(deg$symbol, deg$gene)
rn <- rownames(d$counts); g_sym <- ifelse(rn %in% names(sym) & !is.na(sym[rn]) & sym[rn] != "", sym[rn], rn)
GT <- fread(gtype, header = FALSE, col.names = c("gene_id", "gene_type"))
btype <- GT$gene_type[match(rn, GT$gene_id)]
cat("biotype matched:", mean(!is.na(btype)), "\n")
# 1. technical metrics: done in job 67 (out_genotype_scan_null/technical_*.tsv); not repeated here.
# 2. permutation null
keep <- filterByExpr(d, group = D$dataset); dk <- calcNormFactors(d[keep, , keep.lib.sizes = FALSE])
Y <- cpm(dk, log = TRUE, prior.count = 2); gs <- g_sym[keep]; bt <- btype[keep]; ave_e <- rowMeans(Y)
DROP <- c("FCAMR", "TRIM31", "UBD", "CXCL10", "IL32", "CCL20", "HSD17B13", "PNPLA3")
ok <- !(gs %in% DROP); Y <- Y[ok, ]; gs <- gs[ok]; bt <- bt[ok]; ave_e <- ave_e[ok]
gl <- strsplit(readLines(GMT), "\t"); all_sets <- setNames(lapply(gl, function(x) unique(x[-(1:2)])), sapply(gl, `[`, 1))
HALLG <- unique(unlist(all_sets[grep("^HALLMARK_", names(all_sets))]))
g53 <- fread(file.path(RUN, "out_chemokine_hepatocyte_state/gene_lfc_matched.tsv.gz"))[donors == "masld"][, donors := NULL]
g56 <- fread(file.path(RUN, "out_atlas_state_replication/gene_lfc_matched.tsv.gz"))
SC <- rbind(g53, g56, use.names = TRUE)[mode == "depth" & is.finite(mean_excess) & anchor == "CCL20"]
x <- SC[, .(k = .N, mn = min(mean_excess), mx = max(mean_excess)), by = gene][k == 4]
SETS <- list(CCL20_STATE_UP = gs %in% setdiff(x[mn > 0.25, gene], DROP), CCL20_STATE_DN = gs %in% setdiff(x[mx < -0.25, gene], DROP),
             HALLMARK_ALL = gs %in% HALLG, PROTEIN_CODING = !is.na(bt) & bt == "protein_coding")
pcg <- SETS$PROTEIN_CODING
brk <- unique(quantile(ave_e, seq(0, 1, length.out = 21))); bin <- cut(ave_e, brk, include.lowest = TRUE, labels = FALSE)
mmean <- function(x, s, bg) {
  w <- tabulate(bin[s], max(bin)); bm <- tapply(x[bg], factor(bin[bg], levels = seq_len(max(bin))), mean)
  sum(w * ifelse(is.na(bm), 0, bm)) / sum(w[!is.na(bm)])
}
stats_of <- function(t, refs) {
  out <- c()
  res <- lapply(refs, function(r) { okr <- is.finite(r); e <- rep(NA_real_, length(t)); e[okr] <- residuals(lm.fit(cbind(1, r[okr]), t[okr])); e })
  for (nm in names(SETS)) {
    s <- SETS[[nm]]; bg <- pcg & !s
    out[paste0(nm, "__vs_rest")] <- mean(t[s]) - mean(t[!s])
    if (nm == "PROTEIN_CODING") next
    out[paste0(nm, "__vs_matched_pc")] <- mean(t[s]) - mmean(t, s, bg)
    for (rf in names(res)) { e <- res[[rf]]; okr <- is.finite(e)
      out[paste0(nm, "__resid_", rf, "_vs_matched_pc")] <- mean(e[s & okr]) - mmean(e, s & okr, bg & okr) }
  }
  out
}
VAR <- c("PNPLA3_148M", "HSD17B13_junction")
scan_t <- function(g, C, RY, rss0) {
  rg <- qr.resid(C, g); sg <- sum(rg^2); b <- as.numeric(crossprod(rg, RY)) / sg
  s2 <- (rss0 - b^2 * sg) / (length(g) - C$rank - 1); b / sqrt(s2 / sg)
}
perm_within <- function(g, grp) { for (co in unique(grp)) { j <- which(grp == co); g[j] <- g[j][sample.int(length(j))] }; g }
P <- list()
run_perm <- function(v, idx, B, label, covf) {
  Dv <- D[idx, ]; C <- qr(model.matrix(covf, data = Dv)); RY <- qr.resid(C, t(Y[, idx])); rss0 <- colSums(RY^2)
  g <- Dv[[v]]; obs <- stats_of(scan_t(g, C, RY, rss0), REFS)
  nul <- t(replicate(B, stats_of(scan_t(perm_within(g, Dv$dataset), C, RY, rss0), REFS)))
  saveRDS(nul, file.path(out, paste0("null_", v, "_", label, ".rds")))
  mu <- colMeans(nul); sdv <- apply(nul, 2, sd)
  pp <- sapply(seq_along(obs), function(k) (1 + sum(abs(nul[, k] - mu[k]) >= abs(obs[k] - mu[k]))) / (1 + B))
  data.table(variant = v, scope = label, n = length(idx), stat = names(obs), obs = obs, null_mean = mu, null_sd = sdv,
             z = (obs - mu) / sdv, p_perm = pp, B = B)
}
# fixed reference t vectors
dis_t <- deg$t[match(rn[keep][ok], deg$gene)]
fidx <- which(!is.na(D$fibrosis_stage) & !is.na(D$sex_final))
Cf <- qr(model.matrix(~ dataset + sex_final, data = D[fidx, ])); RYf <- qr.resid(Cf, t(Y[, fidx]))
fib_t <- scan_t(D$fibrosis_stage[fidx], Cf, RYf, colSums(RYf^2)); rm(RYf)
REFS <- list(disease = dis_t, fibrosis = fib_t)
cat("disease t matched:", mean(is.finite(dis_t)), " cor(disease t, fibrosis t):", cor(dis_t, fib_t, use = "complete.obs"), "\n")
for (v in VAR) {
  idx <- which(!is.na(D[[v]]) & !is.na(D$sex_final))
  P[[paste(v, "pooled")]] <- run_perm(v, idx, 10000, "pooled", ~ dataset + sex_final)
  for (co in unique(D$dataset)) {
    j <- idx[D$dataset[idx] == co]
    if (length(unique(D[[v]][j])) < 2) next
    covf <- if (length(unique(D$sex_final[j])) > 1) ~ sex_final else ~ 1
    P[[paste(v, co)]] <- run_perm(v, j, 2000, co, covf)
  }
}
PP <- rbindlist(P); fwrite(PP, file.path(out, "permutation_set_stats.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(PP[, .(variant, scope, n, stat, obs = round(obs, 3), null_mean = round(null_mean, 3), z = round(z, 2), p_perm = signif(p_perm, 2))])
