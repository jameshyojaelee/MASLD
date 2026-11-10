# SS S3g (F-SS7, PRESPEC_SS section 5; exploratory; forks of DIG1 jobs 41a, 58, 10, 42b, 50 and 73). DIG1 expression
# findings against stage_num in SteatoSITE. Directions are the frozen DIG1 signs; one-sided in that direction unless
# marked two-sided; BH over all 14 tests below. Covariates batch + inferred_sex (+ the section's own covariates).
# Units: staged included libraries with inferred_sex. Scores: z over all included libraries, mean over set genes.
#   s45a-c  TNFSF10 < 0, TNFRSF10B > 0, TNFRSF10A > 0: partial Spearman of gene log CPM (TMM, prior 1) with stage_num,
#           ranks residualized on batch + sex (job 41a); genome-wide table with the job 41a filter (log CPM > 1 in half
#           of included libraries), named genes reported whether or not they pass the filter.
#   s38     INFLAMED_HEP (job 58 set; log2 CPM prior 2) ~ stage_num + IMMUNE + HEPMASS + DUCTULAR + batch + sex, OLS HC3, > 0.
#   s9a,b   partial Spearman GNMT vs biliary < 0, GNMT vs CDKN1A < 0 | batch + sex + hepatocyte proxy (job 10 sets, prior 1).
#   s11a-c  partial Spearman with stage_num | batch + sex: biliary > 0, CDKN1A > 0, GNMT < 0 (job 10 stage table).
#   s30a,b  THRB < 0, SLC16A10 < 0, as s45 (job 41a statistic); + hepatocyte proxy beside (job 42b ADJF).
#   s33     injury axis S (job 50: 291 up / 176 down replicated axis genes minus KRT8, KRT18, SQSTM1; job 41a filter;
#           non-hepatocyte proxies stroma, immune, cholangiocyte, endothelial) vs stage_num, partial Spearman | batch
#           + sex + proxies, two-sided.
#   loci    genes of enzyme-only (job 73 class "release", 46 loci) vs disease loci (18): |moderated t| of stage_num from
#           limma voomWithQualityWeights ~ batch + inferred_sex + stage_num (filterByExpr(group = bin)); Wilcoxon
#           rank-sum, two-sided. Locus-mean |t| beside.
#   s45cis  TNFSF10 log2 CPM (prior 2) ~ dosage of the S0d chosen tag allele that sits on the lead G (ALT-raising)
#           haplotype + batch + sex, OLS HC3, > 0. Testable if >= 100 called staged units, >= 2 genotype classes and
#           HWE p >= 1e-4 (calls: s3f_genotype_table.py). Tag -> stage_num beside (two-sided, expected null).
# --dataset-term (774-unit reproduction only): + dataset in every model, z within dataset.
# Unit: genetic individual (S2 Q3). Cross-sectional.
# Usage: Rscript s3g_expression.R --input <S2 dir> --stage <stage tsv> --geno <s3f_genotype_table dir>
#        --targets <out_s0d_targets dir> --out <new dir> [--dataset-term]
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table)})
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k) { i <- match(k, a); if (is.na(i)) NA_character_ else a[i + 1] }
DSET <- "--dataset-term" %in% a
inp <- opt("--input"); stf <- opt("--stage"); gdir <- opt("--geno"); tdir <- opt("--targets"); out <- opt("--out")
stopifnot(!anyNA(c(inp, stf, tdir, out)), !dir.exists(out)); dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SS <- file.path(ROOT, "Analysis/Hypotheses/SS-20261025T135237Z")
stopifnot(system(paste("cd", file.path(SS, "frozen"), "&& sha256sum -c --quiet SHA256SUMS")) == 0)
AXIS <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z/out_length_confound/genome_A2_cross_cohort.tsv.gz")
stopifnot(sub(" .*", "", system2("sha256sum", AXIS, stdout = TRUE)) == "6cf901e4a481f09dadab507052775751725bb3867714c2504862b3f1b5a4348c")
U <- fread(file.path(inp, "units.tsv"))[toupper(as.character(include)) == "TRUE"]
U[, sample_alias := as.character(sample_alias)]
if (!DSET) U[, dataset := "one_cohort"]
st <- fread(cmd = paste("tr -d '\\r' <", shQuote(stf)), colClasses = "character")
U[, stage_num := suppressWarnings(as.numeric(st$nash_crn_kleiner_fibrosis_stage[match(sample_alias, st$sample_alias)]))]
d <- readRDS(file.path(inp, "dge_included.rds")); stopifnot(identical(colnames(d$counts), U$sample_alias))
deg <- fread(file.path(SS, "frozen/release/deg_results.csv"))
rn <- rownames(d$counts); ens <- sub("\\..*$", "", rn); sym <- deg$symbol[match(rn, deg$gene)]
d <- calcNormFactors(d); lc1 <- cpm(d, log = TRUE, prior.count = 1); lc2 <- cpm(d, log = TRUE, prior.count = 2)
COV <- c(if (DSET) "dataset", "batch", "inferred_sex")
sel <- !is.na(U$stage_num) & !is.na(U$inferred_sex)
Z0 <- function(M) t(apply(M, 1, function(v) ave(v, U$dataset, FUN = function(x) (x - mean(x)) / sd(x))))
score <- function(M, genes) { idx <- which(sym %in% genes); colMeans(Z0(M[idx, , drop = FALSE]), na.rm = TRUE) }
used_cov <- function(z, vars) vars[vapply(vars, function(v) length(unique(z[[v]])) > 1, logical(1))]
pcor <- function(z, x, y, vars) {           # partial Spearman (jobs 10, 41a): ranks residualized on vars
  z <- z[complete.cases(z[, c(x, y, vars), with = FALSE])]; vars <- used_cov(z, vars)
  xx <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[x]]), rank(z[[y]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(n = nrow(z), est = r, stat = r * sqrt(df / (1 - r^2)), df = df) }
hc3 <- function(f, z) {
  v <- all.vars(f); z <- z[complete.cases(z[, v, with = FALSE])]
  rhs <- attr(terms(f), "term.labels"); rhs <- c(rhs[1], used_cov(z, rhs[-1]))
  mod <- lm(reformulate(rhs, v[1]), data = z); X <- model.matrix(mod); e <- residuals(mod); h <- hatvalues(mod)
  XtXi <- solve(crossprod(X)); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi; k <- match(rhs[1], colnames(X))
  data.table(n = nobs(mod), est = coef(mod)[k], se = sqrt(V[k, k]), stat = coef(mod)[k] / sqrt(V[k, k]), df = df.residual(mod)) }
ps <- function(stat, df, dir) if (dir == 0) 2 * pt(-abs(stat), df) else pt(stat, df, lower.tail = dir < 0)
TESTS <- list(); add <- function(id, section, desc, dir, r, dig1 = NA_real_) TESTS[[id]] <<- cbind(id = id, section = section, test = desc,
  direction = c("<0", "two-sided", ">0")[dir + 2], r, p = ps(r$stat, r$df, dir), dig1_value = dig1)
D <- U[sel]
# s45, s30: gene-level partial Spearman (job 41a)
keep <- rowMeans(lc1 > 1) >= 0.5
X <- lc1[, sel, drop = FALSE]; q <- qr(model.matrix(reformulate(if (length(used_cov(D, COV))) used_cov(D, COV) else "1"), data = D))
rs <- qr.resid(q, rank(D$stage_num)); RX <- qr.resid(q, apply(X, 1, rank))
r <- as.numeric(crossprod(RX, rs)) / sqrt(colSums(RX^2) * sum(rs^2)); df <- nrow(D) - q$rank - 1
G <- data.table(ensembl = ens, symbol = sym, passes_41a_filter = keep, n = nrow(D), rho = r, t = r * sqrt(df / (1 - r^2)), df = df)
G[, p := 2 * pt(-abs(t), df)][passes_41a_filter == TRUE, bh := p.adjust(p, "BH")]
fwrite(G, file.path(out, "gene_stage_partial_spearman.tsv.gz"), sep = "\t")
dig <- fread(file.path(SS, "frozen/dig1/gene_disease_state.tsv.gz"))
gt <- function(s) { x <- G[symbol == s]; stopifnot(nrow(x) == 1); x[, .(n, est = rho, stat = t, df)] }
for (k in list(c("s45a", "TNFSF10", -1), c("s45b", "TNFRSF10B", 1), c("s45c", "TNFRSF10A", 1), c("s30a", "THRB", -1), c("s30b", "SLC16A10", -1)))
  add(k[1], sub("[abc]$", "", sub("^s", "§", k[1])), paste(k[2], "vs stage_num"), as.integer(k[3]), gt(k[2]), dig[symbol == k[2], rho_fib][1])
# scores
S58 <- list(INFLAMED_HEP = c("FCAMR", "TRIM31", "UBD", "CXCL10", "IL32", "CCL20"), IMMUNE = c("PTPRC", "CD3E", "CD2", "CD68", "CD163"),
            DUCTULAR = c("KRT7", "KRT19", "EPCAM", "SOX9"), HEPMASS = c("ALB", "APOA1", "CYP2E1", "HP"))
for (s in names(S58)) set(U, j = s, value = score(lc2, S58[[s]]))
S10 <- list(biliary = c("SOX9", "EPCAM", "KRT7", "KRT19", "CD24", "SPP1", "KRT23"), p21 = "CDKN1A", GNMT = "GNMT",
            hep_proxy = c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB"))
for (s in names(S10)) set(U, j = paste0("j10_", s), value = score(lc1, S10[[s]]))
# job 50 axis score
lf <- lc1[keep, ]; sf <- sym[keep]; ok <- !is.na(sf) & !(sf %in% sf[duplicated(sf)]); lf <- lf[ok, ]; sf <- sf[ok]
Zf <- Z0(lf); rownames(Zf) <- sf
ax <- fread(AXIS); SG <- c("KRT8", "KRT18", "SQSTM1")
up <- setdiff(intersect(ax[replicated == "up", gene], sf), SG); dn <- setdiff(intersect(ax[replicated == "down", gene], sf), SG)
PROX <- list(stroma = c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "TAGLN", "ACTA2", "MYL9"), immune = c("PTPRC", "CD74", "LYZ", "CD68", "CD14", "CD3E"),
             chol = c("KRT7", "KRT19", "EPCAM"), endo = c("PECAM1", "VWF", "CLEC4G"))
for (k in names(PROX)) set(U, j = paste0(k, "_proxy"), value = colMeans(Zf[setdiff(intersect(PROX[[k]], sf), c(up, dn)), , drop = FALSE]))
U[, S_axis := colMeans(Zf[up, , drop = FALSE]) - colMeans(Zf[dn, , drop = FALSE])]
D <- U[sel]
add("s38", "§38", "INFLAMED_HEP ~ stage_num + IMMUNE + HEPMASS + DUCTULAR", 1,
    hc3(reformulate(c("stage_num", "IMMUNE", "HEPMASS", "DUCTULAR", COV), "INFLAMED_HEP"), D)[, .(n, est, stat, df)])
add("s9a", "§9", "GNMT vs biliary | hep_proxy", -1, pcor(D, "j10_GNMT", "j10_biliary", c(COV, "j10_hep_proxy")), -0.33)
add("s9b", "§9", "GNMT vs CDKN1A | hep_proxy", -1, pcor(D, "j10_GNMT", "j10_p21", c(COV, "j10_hep_proxy")), -0.21)
add("s11a", "§11", "biliary vs stage_num", 1, pcor(D, "j10_biliary", "stage_num", COV), 0.50)
add("s11b", "§11", "CDKN1A vs stage_num", 1, pcor(D, "j10_p21", "stage_num", COV), 0.32)
add("s11c", "§11", "GNMT vs stage_num", -1, pcor(D, "j10_GNMT", "stage_num", COV), -0.29)
add("s33", "§33", "injury axis S vs stage_num | non-hepatocyte proxies", 0,
    pcor(D, "S_axis", "stage_num", c(COV, "stroma_proxy", "immune_proxy", "chol_proxy", "endo_proxy")))
# enzyme-only vs disease loci
D[, bin := factor(fifelse(stage_num <= 1, "F0F1", fifelse(stage_num == 2, "F2", "F3F4")))]
dk <- d[, sel]; X <- model.matrix(reformulate(c(used_cov(D, COV), "stage_num")), data = D)
dk <- dk[filterByExpr(dk, group = D$bin), , keep.lib.sizes = FALSE]; dk <- calcNormFactors(dk)
fit <- eBayes(lmFit(voomWithQualityWeights(dk, X), X))
LT <- data.table(gene = rownames(dk), ensembl = sub("\\..*$", "", rownames(dk)), t = fit$t[, "stage_num"], p = fit$p.value[, "stage_num"])
fwrite(LT, file.path(out, "limma_stage_num_t.tsv.gz"), sep = "\t")
lo <- fread(file.path(SS, "frozen/dig1/loci.tsv")); cg <- fread(file.path(SS, "frozen/dig1/coloc_genes_lead_z.tsv"))
cg <- merge(cg[, .(gene, ensembl, locus)], lo[, .(locus, cls)], by = "locus")[cls %in% c("release", "disease")]
cg <- merge(cg, LT[, .(ensembl, t)], by = "ensembl")
w <- wilcox.test(abs(cg[cls == "release", t]), abs(cg[cls == "disease", t]))
nr <- sum(cg$cls == "release"); nd <- sum(cg$cls == "disease")
TESTS[["loci"]] <- data.table(id = "loci", section = "new", test = "|stage t| genes of enzyme-only vs disease loci (Wilcoxon)", direction = "two-sided",
  n = nr + nd, est = median(abs(cg[cls == "release", t])) - median(abs(cg[cls == "disease", t])), stat = w$statistic, df = NA_real_, p = w$p.value, dig1_value = NA_real_)
lm_ <- cg[, .(mean_abs_t = mean(abs(t)), n_genes = .N), by = .(locus, cls)]
wl <- wilcox.test(lm_[cls == "release", mean_abs_t], lm_[cls == "disease", mean_abs_t])
fwrite(cg, file.path(out, "loci_gene_t.tsv"), sep = "\t")
# TNFSF10 cis
TG <- fread(file.path(tdir, "s1b_targets_annotated.tsv"))[role == "tnfsf10" & grepl("CHOSEN TAG", note)]
cis <- list(testable = FALSE, reason = "no qualifying tag in S0d")
if (nrow(TG) == 1 && !is.na(gdir)) {
  rid <- paste0("TNFSF10_tag_chr3:", TG$pos); onG <- sub(".*allele on lead G haplotype ([ACGT]).*", "\\1", TG$note)
  W <- fread(file.path(gdir, "genotype_calls_wide.tsv")); QC <- fread(file.path(gdir, "site_qc.tsv"))
  g <- W[[rid]][match(U$sample_alias, W$sample_alias)]
  U[, tag_dos := if (onG == TG$alt) g else 2 - g]
  U[, TNFSF10_l2 := lc2[which(sym == "TNFSF10"), ]]
  D <- U[sel]; nc <- sum(!is.na(D$tag_dos)); hw <- QC[rsid == rid, min_cohort_hwe_p]
  cis <- list(testable = nc >= 100 && length(unique(na.omit(D$tag_dos))) >= 2 && isTRUE(hw >= 1e-4), n_called = nc, hwe_p = hw,
              tag = rid, allele_counted = onG, reason = "")
  if (cis$testable) {
    add("s45cis", "§45", "TNFSF10 log2 CPM ~ dosage of the lead-G-haplotype tag allele", 1,
        hc3(reformulate(c("tag_dos", COV), "TNFSF10_l2"), D)[, .(n, est, stat, df)])
    tr <- hc3(reformulate(c("tag_dos", COV), "stage_num"), D)
    cis$tag_to_stage <- paste0("est ", signif(tr$est, 3), " se ", signif(tr$se, 3), " p ", signif(ps(tr$stat, tr$df, 0), 3), " n ", tr$n)
  } else cis$reason <- "fails testability rule"
}
R <- rbindlist(TESTS, fill = TRUE); R[, bh := p.adjust(p, "BH")]
fwrite(R, file.path(out, "f_ss7_tests.tsv"), sep = "\t")
U[, THRB_z := score(lc1, "THRB")][, SLC16A10_z := score(lc1, "SLC16A10")]
D <- U[sel]
BES <- rbind(
  cbind(analysis = "THRB + hep_proxy (job 42b ADJF)", pcor(D, "THRB_z", "stage_num", c(COV, "j10_hep_proxy"))),
  cbind(analysis = "SLC16A10 + hep_proxy (job 42b ADJF)", pcor(D, "SLC16A10_z", "stage_num", c(COV, "j10_hep_proxy"))),
  cbind(analysis = "S axis, no proxies", pcor(D, "S_axis", "stage_num", COV)),
  cbind(analysis = "INFLAMED_HEP ~ stage_num (no score covariates)", hc3(reformulate(c("stage_num", COV), "INFLAMED_HEP"), D)[, .(n, est, stat, df)]),
  cbind(analysis = "s38 without F4 (codes 1-3)", hc3(reformulate(c("stage_num", "IMMUNE", "HEPMASS", "DUCTULAR", COV), "INFLAMED_HEP"), D[stage_num <= 3])[, .(n, est, stat, df)]),
  data.table(analysis = "loci, locus-mean |t| Wilcoxon", n = nrow(lm_), est = NA_real_, stat = wl$statistic, df = NA_real_), fill = TRUE)
BES[, p_two_sided := 2 * pt(-abs(stat), df)]
fwrite(BES, file.path(out, "beside.tsv"), sep = "\t")
writeLines(c(paste("mode:", if (DSET) "dataset-term reproduction" else "SteatoSITE one cohort"), paste("staged units:", sum(sel)),
             paste("axis genes in bulk: up", length(up), "down", length(dn)), paste("loci genes: release", nr, "disease", nd),
             paste("locus-mean Wilcoxon p:", signif(wl$p.value, 3)), paste("TNFSF10 cis:", paste(names(cis), unlist(cis), sep = "=", collapse = "; "))),
           file.path(out, "run_notes.txt"))
fwrite(U, file.path(out, "unit_scores.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220); print(readLines(file.path(out, "run_notes.txt"))); print(R, digits = 3); print(BES, digits = 3)
