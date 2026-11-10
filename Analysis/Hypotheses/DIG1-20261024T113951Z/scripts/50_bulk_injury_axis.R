# DIG1 job 50 (exploratory, post hoc to jobs 46-48). How much of the bulk MASLD transcriptome lines up
# with the injured-hepatocyte state found inside single livers? The axis comes from cells only: genes
# whose within-donor coupling to KRT8+SQSTM1 is BH < 0.05 with the same sign in both GSE298719 Flex and
# GSE256398 snRNA after gene-length control (job 47 A2; 291 up, 176 down). KRT8, KRT18, SQSTM1 are
# excluded. Bulk data are never used to pick genes.
# Bulk: five-cohort fragment-count DGE (774-unit release), log CPM > 1 in half of units, z within cohort.
# Unit score S = mean z(axis up genes) - mean z(axis down genes); S_up and S_down reported separately.
# Covariates: dataset + sex + non-hepatocyte composition proxies (43a markers for stroma, immune,
# cholangiocyte, endothelial; axis genes dropped from proxies). The 43a hepatocyte proxy (ALB, APOA1,
# SERPINA1, TF, HP, AHSG, FGB) is itself secreted-protein output, which the axis lowers, so it is a
# sensitivity only, together with no proxies.
# Primary family (Holm over 4), fixed before any bulk value of S was computed:
#   P1 S, disease vs control, all units (OLS t).
#   P2 S vs NAS, disease units in GSE130970, GSE135251, GSE162694 (partial Spearman).
#   P3 S vs fibrosis stage, disease units in those three plus GSE213621.
#   P4 S vs ballooning grade, GSE130970 disease units.
# Descriptive decomposition: for canonical bulk DEGs (padj < 0.05, |logFC| > 0.5), refit
# z ~ dataset + sex + proxies + group with and without S (S without gene g for axis genes) and report
# the disease beta kept. Null for "is the axis special": 200 scores of random non-axis genes matched to
# each axis gene on bulk disease t decile x mean log CPM decile, same up/down structure; the share of
# disease beta^2 removed is compared over canonical DEGs outside the axis and outside the draw.
# Usage: Rscript 50_bulk_injury_axis.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
y <- readRDS(file.path(B, "merged_dge.rds")); m <- readRDS(file.path(B, "meta_matched.rds"))
m <- m[match(colnames(y), m$sample_id)]; stopifnot(identical(m$sample_id, colnames(y)))
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE, prior.count = 1)
lc <- lc[rowMeans(lc > 1) >= 0.5, ]
ens <- sub("\\..*$", "", rownames(lc)); sym <- an$gene_name[match(ens, an$ensembl_base)]
keep <- !is.na(sym) & !duplicated(sym) & !(sym %in% sym[duplicated(sym)])
lc <- lc[keep, ]; ens <- ens[keep]; sym <- sym[keep]
Z <- t(apply(lc, 1, function(v) ave(v, m$dataset, FUN = function(a) (a - mean(a)) / sd(a))))
rownames(Z) <- sym
ax <- fread(file.path(RUN, "out_length_confound/genome_A2_cross_cohort.tsv.gz"))
SCORE_GENES <- c("KRT8", "KRT18", "SQSTM1")
up <- setdiff(intersect(ax[replicated == "up", gene], sym), SCORE_GENES)
dn <- setdiff(intersect(ax[replicated == "down", gene], sym), SCORE_GENES)
PROX <- list(hep = c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB"),
             stroma = c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "TAGLN", "ACTA2", "MYL9"),
             immune = c("PTPRC", "CD74", "LYZ", "CD68", "CD14", "CD3E"), chol = c("KRT7", "KRT19", "EPCAM"),
             endo = c("PECAM1", "VWF", "CLEC4G"))
used <- list()
for (k in names(PROX)) { gg <- setdiff(intersect(PROX[[k]], sym), c(up, dn)); used[[k]] <- gg
  m[, (paste0(k, "_proxy")) := colMeans(Z[gg, , drop = FALSE])] }
fwrite(data.table(proxy = names(used), genes = vapply(used, paste, "", collapse = ",")), file.path(out, "proxies_used.tsv"), sep = "\t")
fwrite(data.table(gene = c(up, dn), direction = rep(c("up", "down"), c(length(up), length(dn)))), file.path(out, "axis_genes_in_bulk.tsv"), sep = "\t")
cat("axis genes in bulk: up", length(up), "down", length(dn), "\n")
m[, `:=`(S_up = colMeans(Z[up, ]), S_down = colMeans(Z[dn, ]))]
m[, S := S_up - S_down]
m[, grp := factor(group_binary, levels = c("Control", "Disease"))]
COV <- list(nonhep = "+ stroma_proxy + immune_proxy + chol_proxy + endo_proxy", none = "",
            all = "+ hep_proxy + stroma_proxy + immune_proxy + chol_proxy + endo_proxy")

# unit-level tests
grp_t <- function(v, cov, sel = rep(TRUE, nrow(m))) {
  d <- copy(m[sel]); d[["yv"]] <- v[sel]
  f <- lm(as.formula(paste("yv ~ dataset + sex_final + grp", cov)), data = d); s <- summary(f)$coefficients["grpDisease", ]
  data.table(n = nrow(d), beta = s[1], se = s[2], t = s[3], p = s[4]) }
pspear <- function(v, score, sets, cov) {
  sel <- m$group_binary == "Disease" & m$dataset %in% sets & !is.na(suppressWarnings(as.numeric(m[[score]])))
  d <- m[sel]; s <- as.numeric(d[[score]])
  base <- c(if (length(unique(d$dataset)) > 1) "dataset", if (length(unique(d$sex_final)) > 1) "sex_final")
  D <- model.matrix(as.formula(paste("~", paste(c("1", base), collapse = " + "), cov)), data = d)
  q <- qr(D); rs <- qr.resid(q, rank(s)); rv <- qr.resid(q, rank(v[sel]))
  r <- sum(rs * rv) / sqrt(sum(rs^2) * sum(rv^2)); df <- nrow(D) - q$rank - 1
  data.table(n = nrow(d), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
NAS_SETS <- c("GSE130970", "GSE135251", "GSE162694"); FIB_SETS <- c(NAS_SETS, "GSE213621")
unit_rows <- list()
for (cn in names(COV)) for (sc in c("S", "S_up", "S_down")) {
  v <- m[[sc]]
  unit_rows[[length(unit_rows) + 1]] <- rbind(
    cbind(test = "disease_vs_control", grp_t(v, COV[[cn]])[, .(n, est = beta, se, stat = t, p)]),
    cbind(test = "NAS", pspear(v, "nas_score", NAS_SETS, COV[[cn]])[, .(n, est = rho, se = NA_real_, stat = NA_real_, p)]),
    cbind(test = "fibrosis", pspear(v, "fibrosis_stage", FIB_SETS, COV[[cn]])[, .(n, est = rho, se = NA_real_, stat = NA_real_, p)]),
    cbind(test = "ballooning", pspear(v, "ballooning_grade", "GSE130970", COV[[cn]])[, .(n, est = rho, se = NA_real_, stat = NA_real_, p)]))[
      , `:=`(covariates = cn, score = sc)]
}
U <- rbindlist(unit_rows)
fwrite(U, file.path(out, "unit_score_tests.tsv"), sep = "\t")
prim <- U[covariates == "nonhep" & score == "S"][, id := c("P1", "P2", "P3", "P4")]
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
per_ds <- rbindlist(lapply(unique(m$dataset), function(ds) {
  sel <- m$dataset == ds; d <- m[sel]
  if (length(unique(d$grp)) < 2) return(NULL)
  f <- lm(as.formula(paste(if (length(unique(d$sex_final)) > 1) "S ~ sex_final + grp" else "S ~ grp", COV$nonhep)), data = d); s <- summary(f)$coefficients["grpDisease", ]
  data.table(dataset = ds, n_control = sum(d$grp == "Control"), n_disease = sum(d$grp == "Disease"), beta = s[1], se = s[2], p = s[4]) }))
fwrite(per_ds, file.path(out, "disease_effect_per_dataset.tsv"), sep = "\t")

# per-gene decomposition
gfit <- function(Y, X, j) {  # Y genes x units; returns disease beta and t for column j of X
  q <- qr(X); stopifnot(q$rank == ncol(X)); co <- qr.coef(q, t(Y)); res <- t(Y) - X %*% co
  s2 <- colSums(res^2) / (nrow(X) - q$rank); v <- chol2inv(qr.R(q)); jj <- match(j, q$pivot)
  list(beta = co[j, ], t = co[j, ] / sqrt(s2 * v[jj, jj])) }
X0 <- model.matrix(as.formula(paste("~ dataset + sex_final + grp", COV$nonhep)), data = m); jg <- which(colnames(X0) == "grpDisease")
b0 <- gfit(Z, X0, jg)
de <- fread(file.path(B, "deg_results.csv"))[, ensembl := sub("\\..*$", "", gene)]
canon <- de[padj < 0.05 & abs(logFC) > 0.5, ensembl]
is_canon <- ens %in% canon
cat("canonical DEGs in filtered bulk:", sum(is_canon), "\n")
bS <- gfit(Z, cbind(X0, S = m$S), jg)$beta
axis_g <- which(sym %in% c(up, dn))
for (i in axis_g) {  # leave gene i out of its own score
  s_up <- if (sym[i] %in% up) (m$S_up * length(up) - Z[i, ]) / (length(up) - 1) else m$S_up
  s_dn <- if (sym[i] %in% dn) (m$S_down * length(dn) - Z[i, ]) / (length(dn) - 1) else m$S_down
  bS[i] <- gfit(Z[i, , drop = FALSE], cbind(X0, S = s_up - s_dn), jg)$beta
}
G <- data.table(ensembl = ens, symbol = sym, canonical = is_canon, axis = fifelse(sym %in% up, "up", fifelse(sym %in% dn, "down", "")),
                beta_disease = b0$beta, t_disease = b0$t, beta_disease_given_S = bS)
G[, kept := beta_disease_given_S / beta_disease]
# matched null scores
tb <- cut(rank(b0$t), 10, labels = FALSE); mb <- cut(rank(rowMeans(lc)), 10, labels = FALSE); bin <- paste(tb, mb)
pool <- split(which(G$axis == ""), bin[G$axis == ""])
draw <- function(idx) vapply(idx, function(i) { p <- pool[[bin[i]]]; if (length(p) < 5) p <- unlist(pool[startsWith(names(pool), paste0(tb[i], " "))]); p[sample.int(length(p), 1)] }, 1L)
iu <- which(G$axis == "up"); id <- which(G$axis == "down")
eval_set <- G$canonical & G$axis == ""
share_removed <- function(bs, excl) { e <- eval_set & !(seq_len(nrow(G)) %in% excl); 1 - sum(bs[e]^2) / sum(b0$beta[e]^2) }
obs_share <- share_removed(bS, integer())
null_rows <- vector("list", 200)
for (k in 1:200) {
  ru <- draw(iu); rd <- draw(id); sk <- colMeans(Z[ru, , drop = FALSE]) - colMeans(Z[rd, , drop = FALSE])
  bk <- gfit(Z, cbind(X0, S = sk), jg)$beta
  gk <- grp_t(sk, COV$nonhep)
  null_rows[[k]] <- data.table(draw = k, share_removed = share_removed(bk, c(ru, rd)), score_disease_t = gk$t,
                               n_unique = length(unique(c(ru, rd))))
}
N <- rbindlist(null_rows)
fwrite(N, file.path(out, "null_matched_scores.tsv"), sep = "\t")
summ <- data.table(n_units = nrow(m), n_up = length(up), n_down = length(dn), n_canonical_eval = sum(eval_set),
                   share_beta2_removed_axis = obs_share, null_mean = mean(N$share_removed), null_q95 = quantile(N$share_removed, 0.95),
                   p_null = (1 + sum(N$share_removed >= obs_share)) / 201,
                   axis_score_disease_t = prim[id == "P1", stat], null_score_disease_t_mean = mean(N$score_disease_t),
                   median_kept_canonical_nonaxis = G[eval_set, median(kept)],
                   frac_canonical_nonaxis_sign_kept = G[eval_set, mean(sign(beta_disease_given_S) == sign(beta_disease))])
fwrite(summ, file.path(out, "decomposition_summary.tsv"), sep = "\t")
fwrite(G, file.path(out, "gene_disease_beta_with_without_axis.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(prim, digits = 3); print(U[score != "S" | covariates != "nonhep"], digits = 3); print(per_ds, digits = 3); print(summ, digits = 3)
SHOW <- c("GNMT", "MAT1A", "AKR1B10", "FGF21", "ANGPTL8", "INHBE", "IGF1", "IGFALS", "GHR", "SQLE", "HMGCR", "LDLR", "FASN", "SCD", "ELOVL6",
          "CIDEC", "PLIN2", "SPP1", "TNFRSF12A", "IL32", "KRT8", "KRT18", "SQSTM1", "CFH", "C3", "FGB", "GDF15", "CYP1A2", "CYP2C8", "SOX9", "EPCAM",
          "THRSP", "ME1", "PNPLA3", "TM6SF2", "HSD17B13", "CA12", "FMO1", "TRIM31", "SLCO1A2", "CYP1A1", "NFKBIA", "SOD2", "BIRC3", "RELB")
print(G[symbol %in% SHOW][order(kept)], digits = 3)
cat("most kept canonical non-axis (|t| > 8):\n"); print(G[eval_set & abs(t_disease) > 8][order(-abs(kept))][1:30], digits = 3)
cat("least kept canonical non-axis (|t| > 8):\n"); print(G[eval_set & abs(t_disease) > 8][order(abs(kept))][1:30], digits = 3)
