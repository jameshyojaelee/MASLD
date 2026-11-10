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
# Unit: one bulk sample per participant. No ancestry covariates. Cross-sectional.
# Usage: Rscript 67_genotype_scan_null.R <job63_out_dir> <job66_out_dir> <gene_type_tsv> <new_out_dir>
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
# 1. technical metrics
star_dir <- function(ds, s) file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/results", ds, "alignments/star", s)
tech <- rbindlist(lapply(seq_len(nrow(D)), function(i) {
  s <- D$sample_id[i]; p <- star_dir(D$dataset[i], s)
  L <- readLines(file.path(p, paste0(s, ".Log.final.out")))
  val <- function(k) { x <- grep(k, L, fixed = TRUE, value = TRUE)[1]; as.numeric(sub("%", "", trimws(sub(".*\\|", "", x)))) }
  uniq <- val("Uniquely mapped reads number"); spl <- val("Number of splices: Total")
  R <- read.delim(file.path(p, paste0(s, ".ReadsPerGene.out.tab")), header = FALSE)
  assigned <- colSums(R[-(1:4), 2:4]); k <- which.max(assigned)
  data.table(sample_id = s, uniq_pct = val("Uniquely mapped reads %"), multi_pct = val("% of reads mapped to multiple loci"),
             splices_per_unique = spl / uniq, mismatch_pct = val("Mismatch rate per base, %"),
             assigned_frac = assigned[k] / uniq, strand_col = k)
}))
pc <- !is.na(btype) & btype == "protein_coding"
tech$pc_share <- colSums(d$counts[pc, ]) / colSums(d$counts)
tech$log_lib <- log10(colSums(d$counts))
D <- cbind(D, tech[match(D$sample_id, tech$sample_id), -1])
fwrite(D[, c("sample_id", "dataset", names(tech)[-1], "pc_share", "log_lib")], file.path(out, "technical_metrics.tsv"), sep = "\t")
fit <- function(f, dat) {
  v <- all.vars(f); cc <- dat[complete.cases(dat[, v, drop = FALSE]), , drop = FALSE]
  rhs <- attr(terms(f), "term.labels"); rhs <- rhs[!(rhs %in% c("dataset", "sex_final") & sapply(rhs, function(t) t %in% names(cc) && length(unique(cc[[t]])) < 2))]
  f <- reformulate(rhs, response = v[1])
  mod <- lm(f, data = cc); X <- model.matrix(mod); e <- residuals(mod); h <- hatvalues(mod)
  XtXi <- solve(crossprod(X)); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  b <- coef(mod); se <- sqrt(diag(V))
  data.frame(term = names(b), estimate = b, se = se, p = 2 * pt(-abs(b / se), df.residual(mod)), n = nobs(mod), row.names = NULL)
}
VAR <- c("PNPLA3_148M", "HSD17B13_junction")
MET <- c("uniq_pct", "multi_pct", "splices_per_unique", "mismatch_pct", "assigned_frac", "pc_share", "log_lib")
for (m in MET) D[[paste0(m, "_z")]] <- ave(D[[m]], D$dataset, FUN = function(x) as.numeric(scale(x)))
rows <- list()
for (v in VAR) {
  for (m in MET) for (co in c("pooled", "GSE213621")) {
    dd <- if (co == "pooled") D else D[D$dataset == co, ]
    r <- fit(as.formula(paste0(m, "_z ~ ", v, " + dataset + sex_final")), dd); rows[[length(rows) + 1]] <- cbind(r[r$term == v, ], outcome = m, cohort = co, model = "metric")
  }
  for (co in c("pooled", "GSE213621")) {
    dd <- if (co == "pooled") D else D[D$dataset == co, ]
    r <- fit(as.formula(paste0("INFLAMED_HEP ~ ", v, " + dataset + sex_final + ", paste0(MET, "_z", collapse = " + "))), dd)
    rows[[length(rows) + 1]] <- cbind(r[r$term == v, ], outcome = "INFLAMED_HEP", cohort = co, model = "plus_technical")
  }
}
TM <- do.call(rbind, rows); fwrite(TM, file.path(out, "technical_models.tsv"), sep = "\t")
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
stats_of <- function(t) {
  unlist(lapply(names(SETS), function(nm) {
    s <- SETS[[nm]]; a <- mean(t[s]) - mean(t[!s])
    bg <- pcg & !s; w <- tabulate(bin[s], max(bin)); bm <- tapply(t[bg], factor(bin[bg], levels = seq_len(max(bin))), mean)
    m <- sum(w * ifelse(is.na(bm), 0, bm)) / sum(w[!is.na(bm)])
    setNames(c(a, mean(t[s]) - m), paste0(nm, c("__vs_rest", "__vs_matched_pc")))
  }))
}
scan_t <- function(g, C, RY, rss0) {
  rg <- qr.resid(C, g); sg <- sum(rg^2); b <- as.numeric(crossprod(rg, RY)) / sg
  s2 <- (rss0 - b^2 * sg) / (length(g) - C$rank - 1); b / sqrt(s2 / sg)
}
perm_within <- function(g, grp) { for (co in unique(grp)) { j <- which(grp == co); g[j] <- g[j][sample.int(length(j))] }; g }
P <- list()
run_perm <- function(v, idx, B, label, covf) {
  Dv <- D[idx, ]; C <- qr(model.matrix(covf, data = Dv)); RY <- qr.resid(C, t(Y[, idx])); rss0 <- colSums(RY^2)
  g <- Dv[[v]]; obs <- stats_of(scan_t(g, C, RY, rss0))
  nul <- t(replicate(B, stats_of(scan_t(perm_within(g, Dv$dataset), C, RY, rss0))))
  data.table(variant = v, scope = label, n = length(idx), stat = names(obs), obs = obs, null_mean = colMeans(nul), null_sd = apply(nul, 2, sd),
             z = (obs - colMeans(nul)) / apply(nul, 2, sd), p_perm = (1 + colSums(abs(sweep(nul, 2, colMeans(nul))) >= abs(obs - colMeans(nul)))) / (1 + B),
             set_size = sapply(sub("__.*", "", names(obs)), function(nm) sum(SETS[[nm]])))
}
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
print(TM, digits = 3)
print(PP[, .(variant, scope, n, stat, set_size, obs = round(obs, 3), null_mean = round(null_mean, 3), z = round(z, 2), p_perm = signif(p_perm, 2))])
