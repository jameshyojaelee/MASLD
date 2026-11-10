# CONF1 C5 (PRESPEC_C5_C1b.md): per-hepatocyte GNMT, MAT1A, CYP2C19 in GSE289173 (held out).
# Model, permutation and MDE follow c2_methionine_hepatocyte.R; adds the prespecified dilution guard.
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(trailingOnly = TRUE)
IN <- args[1]; OUT <- args[2]
if (is.na(OUT) || dir.exists(OUT)) stop("output directory missing or already exists")
dir.create(OUT)
SEED <- 20261109; B <- 100000L; MIN_HEP <- 30; DETECT <- 0.80
GENES <- c("GNMT", "MAT1A", "CYP2C19")
DIG1_EFFECT <- c(GNMT = 1.89, MAT1A = 0.82, CYP2C19 = 1.72)
MARKERS <- c("COL1A1", "COL3A1", "DCN", "PECAM1", "STAB2", "PTPRC", "CD163", "LYZ")
wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")

pb <- fread(file.path(IN, "pseudobulk_hep.tsv.gz")); meta <- fread(file.path(IN, "donors.tsv"))
sym <- pb$symbol; uniq <- !(sym %in% sym[duplicated(sym)])
meta <- meta[n_hep >= MIN_HEP]
meta[, disease := fifelse(stage_ord == 0, 0, fifelse(stage_ord >= 2, 1, NA_real_))]
cn <- as.matrix(pb[, meta$donor, with = FALSE]); rownames(cn) <- sym
lc <- cpm(calcNormFactors(DGEList(cn)), log = TRUE, prior.count = 1)
detect <- rowMeans(cn > 0)
wr(meta, "donors_eligible.tsv")

hc3 <- function(Y, X, j = 2) {
  XtXi <- solve(crossprod(X)); C <- XtXi %*% t(X); Bh <- C %*% Y; E <- Y - X %*% Bh
  h <- rowSums((X %*% XtXi) * X); W <- E^2 / (1 - h)^2
  se <- sqrt(as.vector(C[j, , drop = FALSE]^2 %*% W)); df <- nrow(X) - ncol(X)
  data.table(beta = Bh[j, ], se = se, df = df, ci_low = Bh[j, ] - qt(0.975, df) * se, ci_high = Bh[j, ] + qt(0.975, df) * se)
}
design <- function(m, xv, covars = c("female", "log10_n")) {
  m <- copy(m); m[, log10_n := log10(n_hep)]
  Z <- cbind(1, as.matrix(m[, covars, with = FALSE])); Z <- Z[, c(TRUE, apply(Z[, -1, drop = FALSE], 2, sd) > 0), drop = FALSE]
  list(Z = Z, X = cbind(Z[, 1], xv = xv, Z[, -1, drop = FALSE]))
}
fit_perm <- function(m, xv, genes, seed, covars = c("female", "log10_n")) {
  d <- design(m, xv, covars)
  Y <- t(lc[genes, m$donor, drop = FALSE])
  f <- cbind(data.table(gene = genes, n = nrow(m), n0 = sum(xv == 0), n1 = sum(xv == 1)), hc3(Y, d$X))
  MZ <- diag(nrow(d$Z)) - d$Z %*% solve(crossprod(d$Z), t(d$Z))
  set.seed(seed)
  D <- vapply(seq_len(B), function(b) sample(xv), numeric(length(xv)))
  Dr <- MZ %*% D; ss <- colSums(Dr^2); bad <- ss < 1e-8
  N <- crossprod(Dr, Y) / ss; N[bad, ] <- NA_real_
  f[, `:=`(null_sd = apply(N, 2, sd, na.rm = TRUE), n_draws = colSums(is.finite(N)),
           p_perm = vapply(seq_along(genes), function(i) { v <- is.finite(N[, i])
             (1 + sum(abs(N[v, i]) >= abs(beta[i]) - 1e-12)) / (sum(v) + 1) }, 1))]
  f[, mde := 2.8 * null_sd]
  f
}

gate <- data.table(gene = GENES, present = GENES %in% sym[uniq])
gate[, detect_frac := ifelse(present, detect[match(gene, sym)], NA_real_)][, testable := present & detect_frac >= DETECT]
wr(gate, "gene_gate.tsv")
wr(rbindlist(lapply(gate[present == TRUE, gene], function(g)
  data.table(gene = g, donor = meta$donor, group = meta$group, logcpm = lc[g, meta$donor]))), "donor_logcpm.tsv")
tg <- gate[testable == TRUE, gene]

# ---- primary
mp <- meta[!is.na(disease)]
prim <- if (length(tg)) fit_perm(mp, mp$disease, tg, SEED) else data.table()
# dilution guard: all unique detected genes, same model, decile of mean ctrl logCPM
ok <- uniq & detect >= DETECT
d <- design(mp, mp$disease); Yall <- t(lc[ok, mp$donor, drop = FALSE])
ball <- as.vector((solve(crossprod(d$X)) %*% t(d$X) %*% Yall)[2, ])
tab <- data.table(symbol = sym[ok], beta = ball, ctrl = rowMeans(lc[ok, mp[disease == 0, donor], drop = FALSE]))
tab[, decile := cut(ctrl, quantile(ctrl, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
tab[, p20 := quantile(beta, 0.2), by = decile][, frac_below := ecdf(beta)(beta), by = decile]
wr(tab, "all_gene_betas.tsv.gz")
if (nrow(prim)) {
  prim <- merge(prim, tab[, .(gene = symbol, decile, decile_p20 = p20, frac_decile_below = frac_below)], by = "gene")
  prim[, `:=`(holm = p.adjust(p_perm, "holm"), dig1_effect = DIG1_EFFECT[gene], passes_guard = beta < decile_p20,
              median_beta_all = median(tab$beta))]
  prim[, decision := fifelse(beta < 0 & holm < 0.05 & passes_guard, "replicated",
                      fifelse(mde <= dig1_effect, "not_replicated", "uninformative"))]
}
wr(prim, "primary.tsv")

# ---- descriptive
desc <- list()
desc$ordinal <- if (length(tg)) fit_perm(meta, as.numeric(meta$stage_ord), tg, SEED + 1)[, analysis := "ordinal_stage_all"] else NULL
m100 <- mp[n_hep >= 100]
desc$min100 <- if (length(tg) && length(unique(m100$disease)) == 2) fit_perm(m100, m100$disease, tg, SEED + 2)[, analysis := "min_100_hep"] else NULL
mk <- intersect(MARKERS, sym[uniq & detect > 0])
desc$markers <- if (length(mk)) fit_perm(mp, mp$disease, mk, SEED + 3)[, analysis := "nonhep_markers_in_hep_pseudobulk"] else NULL
wr(rbindlist(desc, fill = TRUE), "descriptive.tsv")
print(prim); print(rbindlist(desc, fill = TRUE)[, .(analysis, gene, n, beta, ci_low, ci_high, p_perm, mde)])
cat("median beta all genes", median(tab$beta), "n genes", nrow(tab), "\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
