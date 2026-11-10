# CONF1 C2 (PRESPEC_C2.md): per-hepatocyte GNMT, MAT1A, CYP2C19 in four snRNA cohorts.
# Donor pseudobulks only. Model and normalisation follow DIG1 scripts/01_pseudobulk_dig.R section C.
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
if (is.na(OUT) || dir.exists(OUT)) stop("output directory missing or already exists")
dir.create(OUT)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
HYP <- file.path(ROOT, "Analysis/Hypotheses")
V5 <- file.path(HYP, "V5-20261021T122723Z/outputs")
V3 <- file.path(HYP, "V3-20261008T015239Z/pseudobulk")
NF <- file.path(HYP, "NUCLEUS-FIX-20261107/unique_counts")
IN <- c(
  v5_counts = file.path(V5, "job3/pseudobulk_hep.tsv.gz"),
  v5_donors = file.path(V5, "donor_qc.tsv"),
  v5_cells = file.path(V5, "job3/cells_per_donor_set.tsv"),
  v3_counts = file.path(V3, "pseudobulk_hepatocyte.tsv.gz"),
  v3_donors = file.path(V3, "donor_metadata.tsv"),
  v3_cells = file.path(V3, "nuclei_per_donor_lineage.tsv"),
  atlas_counts = file.path(NF, "pseudobulk/Hepatocytes_pseudobulk.csv"),
  atlas_donors = file.path(NF, "metadata/donor_metadata.tsv"),
  atlas_libraries = file.path(NF, "metadata/library_metadata.tsv"),
  atlas_sex = file.path(HYP, "RV1-20261031T155620Z/out-20261102T130905Z/donors.tsv"),
  atlas_F = file.path(ROOT, "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/metadata/documented_F_source.tsv"),
  annotation = file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"))
stopifnot(all(file.exists(IN)))
SEED <- 20261109; B <- 100000L; MIN_HEP <- 30; DETECT <- 0.80
GENES <- c("GNMT", "MAT1A", "CYP2C19")
DIG1_EFFECT <- c(GNMT = 1.89, MAT1A = 0.82, CYP2C19 = 1.72)
REPL <- c("GSE202379", "GSE244832"); ALL4 <- c("GSE298719", "GSE256398", REPL)
wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")
wr(data.table(input = names(IN), path = IN, sha256 = vapply(IN, function(f)
  system2("sha256sum", shQuote(f), stdout = TRUE) |> sub(pattern = " .*", replacement = ""), "")), "inputs_sha256.tsv")

# ---- cohorts: counts (rows = symbol) and donor table -----------------------------------
# donor table columns: donor, stage_ord (0 healthy, 1 steatosis, 2 steatohepatitis/F2-F3, 3 cirrhosis/F4),
# female, n_hep. disease (primary) = 0 if stage_ord 0, 1 if stage_ord >= 2, NA if 1.
annot <- fread(IN["annotation"], select = c("ensembl_base", "gene_name"))
load_v5 <- function() {
  dq <- fread(IN["v5_donors"])[included == TRUE]
  cells <- fread(IN["v5_cells"])[set == "hep"]
  pb <- fread(IN["v5_counts"])
  meta <- data.table(donor = dq$gsm, female = as.integer(dq$sex_f),
                     stage_ord = ifelse(dq$disease == 0, 0, ifelse(dq$fibrosis >= 4, 3, 2)),
                     n_hep = cells$n_cells[match(dq$gsm, cells$gsm)])
  stopifnot(all(dq$disease == 0 | dq$fibrosis >= 2))
  list(meta = meta, counts = as.matrix(pb[, meta$donor, with = FALSE]),
       symbol = annot$gene_name[match(sub("\\..*$", "", pb$gene_id), annot$ensembl_base)])
}
load_v3 <- function() {
  m <- fread(IN["v3_donors"]); n3 <- fread(IN["v3_cells"]); pb <- fread(IN["v3_counts"])
  meta <- data.table(donor = m$donor, female = as.integer(m$sex == "female"), stage_ord = m$stage_ordinal,
                     n_hep = n3$hepatocyte[match(m$donor, n3$donor)])
  list(meta = meta, counts = as.matrix(pb[, meta$donor, with = FALSE]), symbol = pb$symbol)
}
atlas_pb <- NULL
load_atlas <- function(ds) {
  if (is.null(atlas_pb)) atlas_pb <<- fread(IN["atlas_counts"])
  dm <- fread(IN["atlas_donors"])[dataset == ds & toupper(as.character(exclude)) == "FALSE"]
  lib <- fread(IN["atlas_libraries"])[donor %in% dm$donor]
  sex <- fread(IN["atlas_sex"], select = c("donor", "sex"))
  stopifnot(all(dm$donor %in% sex$donor), all(lib$sample %in% names(atlas_pb)),
            all(sex$sex[match(dm$donor, sex$donor)] %in% c("Male", "Female")))
  counts <- vapply(dm$donor, function(d) rowSums(as.matrix(atlas_pb[, lib[donor == d, sample], with = FALSE])),
                   numeric(nrow(atlas_pb)))
  meta <- data.table(donor = dm$donor, female = as.integer(sex$sex[match(dm$donor, sex$donor)] == "Female"),
                     stage_ord = dm$stage_ordinal,
                     n_hep = vapply(dm$donor, function(d) sum(lib[donor == d][["n_Hepatocytes"]]), numeric(1)))
  list(meta = meta, counts = counts, symbol = atlas_pb$gene)
}
prep <- function(x, cohort) {
  x$meta[, cohort := cohort]
  x$meta[, disease := fifelse(stage_ord == 0, 0, fifelse(stage_ord >= 2, 1, NA_real_))]
  stopifnot(!anyNA(x$meta[, .(female, stage_ord, n_hep)]))
  x$meta[, eligible := n_hep >= MIN_HEP]
  keep <- x$meta$eligible
  cn <- x$counts[, keep, drop = FALSE]
  x$lc <- cpm(calcNormFactors(DGEList(cn)), log = TRUE, prior.count = 1)
  x$detect <- rowMeans(cn > 0)
  x$unique <- !is.na(x$symbol) & !(x$symbol %in% x$symbol[duplicated(x$symbol)])
  colnames(x$lc) <- x$meta$donor[keep]
  x
}
CO <- list(GSE298719 = prep(load_v5(), "GSE298719"), GSE256398 = prep(load_v3(), "GSE256398"),
           GSE202379 = prep(load_atlas("GSE202379"), "GSE202379"), GSE244832 = prep(load_atlas("GSE244832"), "GSE244832"))
rm(atlas_pb)
wr(rbindlist(lapply(CO, `[[`, "meta"), use.names = TRUE), "donors.tsv")
row_of <- function(x, g) { i <- which(x$symbol == g & x$unique); if (length(i) == 1) i else NA_integer_ }
gate <- rbindlist(lapply(names(CO), function(k) rbindlist(lapply(GENES, function(g) {
  i <- row_of(CO[[k]], g)
  data.table(cohort = k, gene = g, present = !is.na(i), detect_frac = if (is.na(i)) NA_real_ else CO[[k]]$detect[i])
}))))
gate[, testable := present & detect_frac >= DETECT]
wr(gate, "gene_gate.tsv")
wr(rbindlist(lapply(names(CO), function(k) { x <- CO[[k]]
  rbindlist(lapply(gate[cohort == k & present == TRUE, gene], function(g)
    data.table(cohort = k, gene = g, donor = colnames(x$lc), logcpm = x$lc[row_of(x, g), ])))
})), "donor_logcpm.tsv")

# ---- model ----------------------------------------------------------------------------
hc3 <- function(Y, X, j = 2) {
  XtXi <- solve(crossprod(X)); C <- XtXi %*% t(X); Bh <- C %*% Y; E <- Y - X %*% Bh
  h <- rowSums((X %*% XtXi) * X); W <- E^2 / (1 - h)^2
  se <- sqrt(as.vector(C[j, , drop = FALSE]^2 %*% W)); df <- nrow(X) - ncol(X)
  data.table(beta = Bh[j, ], se = se, df = df, ci_low = Bh[j, ] - qt(0.975, df) * se, ci_high = Bh[j, ] + qt(0.975, df) * se)
}
# One analysis = a choice of predictor, donors and covariates in each cohort.
# spec$x(meta) returns the predictor (NA drops the donor); spec$covars names covariate columns.
cohort_fit <- function(k, spec, genes, seed_offset) {
  x <- CO[[k]]
  m <- x$meta[eligible == TRUE & n_hep >= spec$min_hep]
  m[, xv := spec$x(m)]; m <- m[!is.na(xv)]
  m[, log10_n := log10(n_hep)]
  Z <- cbind(1, as.matrix(m[, spec$covars, with = FALSE]))
  Z <- Z[, c(TRUE, apply(Z[, -1, drop = FALSE], 2, sd) > 0), drop = FALSE]   # drop a constant covariate
  rows <- vapply(genes, function(g) if (gate[cohort == k & gene == g, testable]) row_of(x, g) else NA_integer_, 1L)
  ok <- !is.na(rows)
  out <- list(cohort = k, n = nrow(m), n0 = sum(m$xv == 0), n1 = sum(m$xv == 1), w = sum((m$xv - mean(m$xv))^2),
              genes = genes[ok], fit = NULL, null = NULL, n_dropped = 0L)
  if (!any(ok)) return(out)
  Y <- t(x$lc[rows[ok], m$donor, drop = FALSE]); colnames(Y) <- genes[ok]
  out$fit <- cbind(data.table(cohort = k, gene = genes[ok], n = nrow(m)), hc3(Y, cbind(Z[, 1], xv = m$xv, Z[, -1, drop = FALSE])))
  MZ <- diag(nrow(Z)) - Z %*% solve(crossprod(Z), t(Z))
  set.seed(SEED + seed_offset)
  D <- vapply(seq_len(B), function(b) sample(m$xv), numeric(nrow(m)))
  Dr <- MZ %*% D; ss <- colSums(Dr^2); bad <- ss < 1e-8
  N <- crossprod(Dr, Y) / ss; N[bad, ] <- NA_real_
  out$null <- N; out$n_dropped <- sum(bad)
  out
}
combine <- function(fits, set, g) {
  f <- Filter(function(z) z$cohort %in% set && g %in% z$genes, fits)
  if (!length(f)) return(NULL)
  w <- vapply(f, `[[`, 1, "w"); b <- vapply(f, function(z) z$fit[gene == g, beta], 1)
  Tn <- Reduce(`+`, Map(function(z, wi) wi * z$null[, g], f, w)) / sum(w)
  Tobs <- sum(w * b) / sum(w); v <- is.finite(Tn)
  data.table(gene = g, cohorts = paste(vapply(f, `[[`, "", "cohort"), collapse = ";"), n_cohorts = length(f),
             n_donors = sum(vapply(f, `[[`, 1, "n")), T = Tobs, n_negative = sum(b < 0),
             null_sd = sd(Tn[v]), mde = 2.8 * sd(Tn[v]), p_perm = (1 + sum(abs(Tn[v]) >= abs(Tobs) - 1e-12)) / (sum(v) + 1),
             n_draws = sum(v))
}
run <- function(spec) lapply(seq_along(ALL4), function(i) cohort_fit(ALL4[i], spec, GENES, i))
cohort_p <- function(fits) rbindlist(lapply(fits, function(z) if (is.null(z$fit)) NULL else
  cbind(z$fit, n0 = z$n0, n1 = z$n1, weight = z$w, n_draws_dropped = z$n_dropped,
        p_perm = vapply(z$genes, function(g) { nl <- z$null[, g]; v <- is.finite(nl)
          (1 + sum(abs(nl[v]) >= abs(z$fit[gene == g, beta]) - 1e-12)) / (sum(v) + 1) }, 1))))

COV <- c("female", "log10_n")
PRIMARY <- list(x = function(m) m$disease, covars = COV, min_hep = MIN_HEP)
fp <- run(PRIMARY)
pc <- cohort_p(fp)
wr(pc, "per_cohort.tsv")

# ---- primary: replication set ----------------------------------------------------------
prim <- rbindlist(lapply(GENES, function(g) combine(fp, REPL, g)))
prim[, `:=`(analysis = "primary_replication", holm = p.adjust(p_perm, "holm"), dig1_effect = DIG1_EFFECT[gene])]
prim[, decision := fifelse(T < 0 & holm < 0.05 & n_negative == n_cohorts & n_cohorts == length(REPL), "replicated",
                    fifelse(mde <= dig1_effect, "not_replicated", "uninformative"))]
four <- rbindlist(lapply(GENES, function(g) combine(fp, ALL4, g)))[, analysis := "four_cohort_descriptive"]
disc <- rbindlist(lapply(GENES, function(g) combine(fp, setdiff(ALL4, REPL), g)))[, analysis := "dig1_cohorts_descriptive"]
wr(rbindlist(list(prim, four, disc), fill = TRUE), "combined.tsv")

# ---- heterogeneity (inverse-variance with HC3; descriptive) ------------------------------
wr(pc[, { wi <- 1 / se^2; fe <- sum(wi * beta) / sum(wi); Q <- sum(wi * (beta - fe)^2); k <- .N
          tau2 <- max(0, (Q - (k - 1)) / (sum(wi) - sum(wi^2) / sum(wi))); wr_ <- 1 / (se^2 + tau2)
          list(k = k, fixed = fe, fixed_se = sqrt(1 / sum(wi)), Q = Q, Q_p = if (k > 1) pchisq(Q, k - 1, lower.tail = FALSE) else NA_real_,
               I2 = if (k > 1) max(0, (Q - (k - 1)) / Q) else NA_real_, tau2 = tau2,
               random = sum(wr_ * beta) / sum(wr_), random_se = sqrt(1 / sum(wr_))) }, by = gene], "heterogeneity.tsv")
wr(rbindlist(lapply(ALL4, function(drop) rbindlist(lapply(GENES, function(g) combine(fp, setdiff(ALL4, drop), g)))[, left_out := drop])),
   "leave_one_out.tsv")

# ---- sensitivity and ordinal slope -------------------------------------------------------
SENS <- list(
  ordinal_stage_all_donors = list(x = function(m) as.numeric(m$stage_ord), covars = COV, min_hep = MIN_HEP),
  no_nuclei_covariate = list(x = function(m) m$disease, covars = "female", min_hep = MIN_HEP),
  min_50_nuclei = list(x = function(m) m$disease, covars = COV, min_hep = 50),
  steatosis_as_disease = list(x = function(m) as.numeric(m$stage_ord >= 1), covars = COV, min_hep = MIN_HEP))
wr(rbindlist(lapply(names(SENS), function(s) { f <- run(SENS[[s]])
  rbindlist(list(
    rbindlist(lapply(GENES, function(g) combine(f, REPL, g)))[, set := "replication"],
    rbindlist(lapply(GENES, function(g) combine(f, ALL4, g)))[, set := "four_cohort"],
    cohort_p(f)[, .(gene, cohorts = cohort, n_cohorts = 1L, n_donors = n, T = beta, p_perm, set = "single_cohort")]),
    fill = TRUE)[, analysis := s] })), "sensitivity.tsv")

# ---- GSE202379: slope per documented fibrosis stage --------------------------------------
Fs <- fread(IN["atlas_F"], select = c("donor", "F_stage_documented"))
x <- CO$GSE202379
x$meta[, F_doc := suppressWarnings(as.numeric(Fs$F_stage_documented[match(donor, Fs$donor)]))]
CO$GSE202379 <- x
fF <- cohort_fit("GSE202379", list(x = function(m) m$F_doc, covars = COV, min_hep = MIN_HEP), GENES, 99)
wr(cohort_p(list(fF)), "f_stage_GSE202379.tsv")

# ---- scale check: all genes, primary contrast --------------------------------------------
scale_set <- function(set, targets) {
  per <- lapply(set, function(k) { x <- CO[[k]]
    m <- x$meta[eligible == TRUE & !is.na(disease)]; m[, log10_n := log10(n_hep)]
    Z <- cbind(1, m$female, m$log10_n); Z <- Z[, c(TRUE, apply(Z[, -1, drop = FALSE], 2, sd) > 0), drop = FALSE]
    dr <- (diag(nrow(Z)) - Z %*% solve(crossprod(Z), t(Z))) %*% m$disease
    ok <- x$unique & x$detect >= DETECT
    Y <- x$lc[ok, m$donor, drop = FALSE]
    data.table(symbol = x$symbol[ok], beta = as.vector(Y %*% dr) / sum(dr^2),
               healthy = rowMeans(Y[, m$disease == 0, drop = FALSE]), w = sum((m$disease - mean(m$disease))^2))
  })
  common <- Reduce(intersect, lapply(per, `[[`, "symbol"))
  tab <- data.table(symbol = common,
    T = Reduce(`+`, lapply(per, function(p) p$w[1] * p$beta[match(common, p$symbol)])) / sum(vapply(per, function(p) p$w[1], 1)),
    healthy = Reduce(`+`, lapply(per, function(p) p$healthy[match(common, p$symbol)])) / length(per))
  tab[, decile := cut(healthy, quantile(healthy, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
  tab[symbol %in% targets, .(gene = symbol, set = paste(set, collapse = ";"), n_genes_common = nrow(tab), T, healthy_logcpm = healthy,
    decile, n_in_decile = vapply(decile, function(d) sum(tab$decile == d), 1L),
    frac_decile_below = mapply(function(t, d) mean(tab$T[tab$decile == d] < t), T, decile),
    frac_all_below = vapply(T, function(t) mean(tab$T < t), 1), median_T_all = median(tab$T))]
}
testable_in <- function(set) GENES[vapply(GENES, function(g) all(gate[gene == g & cohort %in% set, testable]), TRUE)]
wr(rbindlist(list(scale_set(REPL, testable_in(REPL)), scale_set(ALL4, testable_in(ALL4)),
                  scale_set(setdiff(ALL4, "GSE298719"), testable_in(setdiff(ALL4, "GSE298719"))))), "scale_check.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
