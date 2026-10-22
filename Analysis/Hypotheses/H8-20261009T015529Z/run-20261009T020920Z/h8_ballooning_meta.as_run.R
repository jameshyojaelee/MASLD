#!/usr/bin/env Rscript
# H8: transcripts that track hepatocyte ballooning grade after adjustment for
# steatosis grade, fibrosis stage, lobular inflammation grade and sex, fitted
# per cohort with limma-voom (voomWithQualityWeights) and combined across three
# cohorts by DerSimonian-Laird random-effects meta-analysis.
# Prespecification: Analysis/Hypotheses/H8-20261009T015529Z/PRESPEC_H8.md
#
# Usage: Rscript h8_ballooning_meta.R <new_out_dir> <n_perm> <n_cores>
# The output directory must not exist. Seed 20260923.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(jsonlite)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
out <- args[[1]]
n_perm <- as.integer(args[[2]])
n_cores <- as.integer(args[[3]])
if (dir.exists(out)) stop("Output directory exists: ", out)
dir.create(out, recursive = TRUE)
SEED <- 20260923L
t_start <- Sys.time()

repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
kam <- file.path(repo, "RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T020627Z/external/kamzolas_v1/data")
inputs <- list(
  gse130970_counts = list(
    path = file.path(repo, "figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds"),
    sha = "f931e340546d679dc69dd0d6f712f4e73b1d805e80c741b22d3ed25498c68675"),
  gse130970_meta = list(
    path = file.path(repo, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds"),
    sha = "429e27376ad421f5f5fae19c75f5e760030e727fab65c5768ef562a9939edfc3"),
  gse267145_counts = list(
    path = file.path(repo, "Analysis/MASLD_Model_Benchmark/executions/gse267145-matrix-audit-21064333/raw/GSE269412_RNA.txt.gz"),
    sha = "fed9333cb63ecca817c8dbc6f318c8b0795852df2567685bf509b90e3e1690f3"),
  gse267145_meta = list(
    path = file.path(repo, "Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"),
    sha = "60bc52b86d7c7547c1371927513ef848c9062a7d3b6e17c032c04862d9d30b91"),
  ucam_counts = list(
    path = file.path(kam, "merged_counts.csv"),
    sha = "7e7b4a4616565d00d06ec42cf4849b47fe321af264b0bfc580e0f159f50894cc"),
  ucam_meta = list(
    path = file.path(kam, "metadata.csv"),
    sha = "4d940394042a5f74fb0246861d096f0290e2777b9460fc6f942f97cf7e79dbb9"),
  gene_metadata = list(
    path = file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
    sha = NA_character_)
)
sha256 <- function(p) sub(" .*", "", system2("sha256sum", shQuote(p), stdout = TRUE))
input_log <- rbindlist(lapply(names(inputs), function(k) {
  obs <- sha256(inputs[[k]]$path)
  if (!is.na(inputs[[k]]$sha) && obs != inputs[[k]]$sha) stop("Checksum mismatch: ", k)
  data.table(input = k, path = inputs[[k]]$path, sha256 = obs)
}))
fwrite(input_log, file.path(out, "inputs_sha256.tsv"), sep = "\t")

gene_meta <- fread(inputs$gene_metadata$path, select = c("gene_name", "ensembl_base"))
gene_meta <- unique(gene_meta[!grepl("_PAR_Y", ensembl_base)], by = "ensembl_base")

strip_version <- function(m) {
  ids <- rownames(m)
  m <- m[!grepl("_PAR_Y", ids), , drop = FALSE]
  base <- sub("\\..*$", "", rownames(m))
  dup <- base[duplicated(base)]
  keep <- !(base %in% dup)
  m <- m[keep, , drop = FALSE]
  rownames(m) <- base[keep]
  m
}

# ------------------------------------------------------------------ cohorts
cohorts <- list()

meta130 <- as.data.table(readRDS(inputs$gse130970_meta$path))[dataset == "GSE130970"]
counts130 <- readRDS(inputs$gse130970_counts$path)
stopifnot(all(meta130$sample_id %in% colnames(counts130)))
cohorts$GSE130970 <- list(
  counts = strip_version(counts130[, meta130$sample_id, drop = FALSE]),
  meta = meta130[, .(person = sample_id, ballooning = as.integer(ballooning_grade),
                     steatosis = as.integer(steatosis_grade), fibrosis = as.integer(fibrosis_stage),
                     inflammation = as.integer(lobular_grade), sex = inferred_sex)])
rm(counts130); invisible(gc())

pj <- fread(inputs$gse267145_meta$path)
x269 <- fread(inputs$gse267145_counts$path)
m269 <- as.matrix(x269[, -1]); rownames(m269) <- x269[[1]]; rm(x269)
stopifnot(all(pj$participant_id %in% colnames(m269)), !anyDuplicated(pj$participant_id))
cohorts$GSE267145 <- list(
  counts = strip_version(m269[, pj$participant_id, drop = FALSE]),
  meta = pj[, .(person = participant_id, ballooning = as.integer(ballooning),
                steatosis = as.integer(steatosis), fibrosis = as.integer(fibrosis),
                inflammation = as.integer(lobular_inflammation), sex = sex)])
rm(m269)

um <- fread(inputs$ucam_meta$path)[Dataset == "UCAM"]
uc <- fread(inputs$ucam_counts$path)
mu <- as.matrix(uc[, -1]); rownames(mu) <- uc[[1]]; rm(uc)
stopifnot(all(um[["Sample name"]] %in% colnames(mu)), !anyDuplicated(um[["Sample name"]]))
cohorts[["E-MTAB-9815"]] <- list(
  counts = strip_version(mu[, um[["Sample name"]], drop = FALSE]),
  meta = um[, .(person = `Sample name`, ballooning = as.integer(Ballooning),
                steatosis = as.integer(Steatosis), fibrosis = as.integer(Fibrosis),
                inflammation = as.integer(Inflammation), sex = Sex)])
rm(mu); invisible(gc())

cohort_order <- c("GSE130970", "GSE267145", "E-MTAB-9815")
for (k in cohort_order) {
  m <- cohorts[[k]]$meta
  ok <- complete.cases(m[, .(ballooning, steatosis, fibrosis, inflammation, sex)])
  if (!all(ok)) message(k, ": dropping ", sum(!ok), " persons with missing covariates")
  cohorts[[k]]$meta <- m[ok]
  cohorts[[k]]$counts <- cohorts[[k]]$counts[, m$person[ok], drop = FALSE]
  stopifnot(!anyDuplicated(cohorts[[k]]$meta$person))
  cohorts[[k]]$meta[, sex := factor(sex, levels = c("F", "M"))]
  stopifnot(!anyNA(cohorts[[k]]$meta$sex))
}

# ------------------------------------------------------ qualification table
qual <- rbindlist(lapply(cohort_order, function(k) {
  m <- cohorts[[k]]$meta
  data.table(cohort = k, n_persons = nrow(m), n_genes_input = nrow(cohorts[[k]]$counts),
             balloon_0 = sum(m$ballooning == 0), balloon_1 = sum(m$ballooning == 1),
             balloon_2 = sum(m$ballooning == 2), n_balloon_ge1 = sum(m$ballooning >= 1),
             n_female = sum(m$sex == "F"), n_male = sum(m$sex == "M"),
             qualifies = nrow(m) >= 40 & sum(m$ballooning >= 1) >= 10)
}))
fwrite(qual, file.path(out, "cohort_qualification.tsv"), sep = "\t")
if (sum(qual$qualifies) < 3) {
  writeLines("Fewer than three cohorts qualify; stop.", file.path(out, "STOPPED.txt"))
  quit(save = "no", status = 0)
}

# ------------------------------------------------------- multicollinearity
colin <- rbindlist(lapply(cohort_order, function(k) {
  m <- cohorts[[k]]$meta
  r2a <- summary(lm(ballooning ~ steatosis + inflammation + fibrosis, data = m))$r.squared
  r2b <- summary(lm(ballooning ~ steatosis + inflammation + fibrosis + sex, data = m))$r.squared
  sp <- function(v) suppressWarnings(cor(m$ballooning, m[[v]], method = "spearman"))
  data.table(cohort = k, r2_vs_steat_infl_fib = r2a, vif = 1 / (1 - r2a),
             independent_fraction = 1 - r2a, r2_with_sex = r2b, independent_fraction_with_sex = 1 - r2b,
             spearman_steatosis = sp("steatosis"), spearman_inflammation = sp("inflammation"),
             spearman_fibrosis = sp("fibrosis"),
             little_power_flag = (1 - r2a) < 0.30)
}))
fwrite(colin, file.path(out, "collinearity.tsv"), sep = "\t")

# ---------------------------------------------------------------- fitting
norm_dge <- function(counts, keep) {
  dge <- DGEList(counts[keep, , drop = FALSE])
  dge <- calcNormFactors(dge, method = "RLE")
  method <- "RLE"
  if (any(!is.finite(dge$samples$norm.factors))) {
    dge <- calcNormFactors(dge, method = "TMM"); method <- "TMM"
  }
  list(dge = dge, method = method)
}

fit_terms <- function(dge, design, terms) {
  v <- voomWithQualityWeights(dge, design)
  fit <- eBayes(lmFit(v, design))
  lapply(setNames(terms, terms), function(tm) {
    data.table(gene_id = rownames(fit$coefficients), logFC = unname(fit$coefficients[, tm]),
               se = unname(fit$stdev.unscaled[, tm] * sqrt(fit$s2.post)),
               t = unname(fit$t[, tm]), p = unname(fit$p.value[, tm]), AveExpr = unname(fit$Amean))
  })
}

design_primary <- function(m) model.matrix(~ ballooning + steatosis + fibrosis + inflammation + sex, data = m)

merge_sparse <- function(x, min_n = 3L) {
  x <- as.integer(x)
  repeat {
    tab <- table(x); lv <- as.integer(names(tab))
    if (length(lv) <= 1L || all(tab >= min_n)) break
    i <- which(tab < min_n)[1]
    target <- if (i > 1L) lv[i - 1L] else lv[i + 1L]
    x[x == lv[i]] <- target
  }
  factor(x)
}
design_categorical <- function(m) {
  d <- model.matrix(~ ballooning + merge_sparse(steatosis) + merge_sparse(fibrosis) +
                      merge_sparse(inflammation) + sex, data = m)
  ne <- nonEstimable(d)
  if (!is.null(ne)) d <- d[, setdiff(colnames(d), ne), drop = FALSE]
  d
}
design_unadjusted <- function(m) model.matrix(~ ballooning + sex, data = m)

prep <- list()
for (k in cohort_order) {
  m <- cohorts[[k]]$meta
  d <- design_primary(m)
  keep <- filterByExpr(DGEList(cohorts[[k]]$counts), design = d)
  nd <- norm_dge(cohorts[[k]]$counts, keep)
  prep[[k]] <- list(dge = nd$dge, norm_method = nd$method, meta = m)
  cohorts[[k]]$counts <- NULL
}
invisible(gc())

terms <- c("ballooning", "fibrosis", "steatosis")
obs <- list()
for (k in cohort_order) {
  obs[[k]] <- fit_terms(prep[[k]]$dge, design_primary(prep[[k]]$meta), terms)
  for (tm in terms) {
    tb <- obs[[k]][[tm]]
    tb[, q_within := p.adjust(p, "BH")]
    tb[, gene_name := gene_meta$gene_name[match(gene_id, gene_meta$ensembl_base)]]
    obs[[k]][[tm]] <- tb
    fwrite(tb, file.path(out, sprintf("per_cohort_%s_%s.tsv.gz", k, tm)), sep = "\t")
  }
}

meta_genes <- Reduce(intersect, lapply(cohort_order, function(k) rownames(prep[[k]]$dge)))
meta_genes <- sort(meta_genes)
if (length(meta_genes) < 5000L) {
  writeLines(sprintf("Meta family has %d genes (< 5000); stop.", length(meta_genes)), file.path(out, "STOPPED.txt"))
  quit(save = "no", status = 0)
}

# Vectorized DerSimonian-Laird random-effects meta-analysis.
dl_meta <- function(Y, V) {
  W <- 1 / V
  sw <- rowSums(W)
  ybar <- rowSums(W * Y) / sw
  Q <- rowSums(W * (Y - ybar)^2)
  dfq <- ncol(Y) - 1
  C <- sw - rowSums(W^2) / sw
  tau2 <- pmax(0, (Q - dfq) / C)
  Ws <- 1 / (V + tau2)
  mu <- rowSums(Ws * Y) / rowSums(Ws)
  se <- sqrt(1 / rowSums(Ws))
  z <- mu / se
  p <- 2 * pnorm(-abs(z))
  data.table(mu = mu, se_mu = se, z = z, p_meta = p, q_meta = p.adjust(p, "BH"),
             tau2 = tau2, Q = Q, I2 = pmax(0, (Q - dfq) / Q))
}

assemble <- function(tabs) {
  # tabs: named list (cohort) of data.tables with gene_id, logFC, se, p
  Y <- sapply(cohort_order, function(k) tabs[[k]]$logFC[match(meta_genes, tabs[[k]]$gene_id)])
  S <- sapply(cohort_order, function(k) tabs[[k]]$se[match(meta_genes, tabs[[k]]$gene_id)])
  P <- sapply(cohort_order, function(k) tabs[[k]]$p[match(meta_genes, tabs[[k]]$gene_id)])
  M <- dl_meta(Y, S^2)
  same_sign <- rowSums(Y > 0) == ncol(Y) | rowSums(Y < 0) == ncol(Y)
  n_nominal <- rowSums(P < 0.05)
  M[, `:=`(same_sign = same_sign, n_nominal = n_nominal,
           replicated = q_meta < 0.05 & same_sign & n_nominal >= 2L)]
  list(Y = Y, S = S, P = P, M = M)
}

write_meta <- function(a, file) {
  dt <- data.table(gene_id = meta_genes,
                   gene_name = gene_meta$gene_name[match(meta_genes, gene_meta$ensembl_base)])
  for (j in seq_along(cohort_order)) {
    k <- cohort_order[j]
    dt[, (paste0("logFC_", k)) := a$Y[, j]]
    dt[, (paste0("se_", k)) := a$S[, j]]
    dt[, (paste0("p_", k)) := a$P[, j]]
  }
  dt <- cbind(dt, a$M)
  setorder(dt, p_meta)
  fwrite(dt, file, sep = "\t")
  dt
}

meta_obs <- list()
for (tm in terms) {
  a <- assemble(lapply(obs, `[[`, tm))
  meta_obs[[tm]] <- write_meta(a, file.path(out, sprintf("meta_%s.tsv.gz", tm)))
}
rep_b <- meta_obs$ballooning[replicated == TRUE]
fwrite(rep_b, file.path(out, "replicated_ballooning_genes.tsv"), sep = "\t")
N_obs <- nrow(rep_b)

# --------------------------------------------------------- positive controls
pc_genes <- data.table(term = c("fibrosis", "fibrosis", "steatosis", "steatosis"),
                       gene_id = c("ENSG00000108821", "ENSG00000164692", "ENSG00000147872", "ENSG00000187288"),
                       symbol = c("COL1A1", "COL1A2", "PLIN2", "CIDEC"))
pc <- rbindlist(lapply(seq_len(nrow(pc_genes)), function(i) {
  r <- meta_obs[[pc_genes$term[i]]][gene_id == pc_genes$gene_id[i]]
  if (nrow(r) == 0L) return(data.table(pc_genes[i], in_meta_family = FALSE))
  cbind(pc_genes[i], in_meta_family = TRUE,
        r[, c(grep("^logFC_|^p_", names(r), value = TRUE), "mu", "se_mu", "q_meta", "replicated"), with = FALSE])
}), fill = TRUE)
pc[, gene_pass := in_meta_family & !is.na(mu) & mu > 0 & q_meta < 0.05]
fib_pass <- all(pc[term == "fibrosis", gene_pass])
ste_pass <- any(pc[term == "steatosis", gene_pass])
fwrite(pc, file.path(out, "positive_controls.tsv"), sep = "\t")

# ------------------------------------------------------- secondary analyses
secondary <- list()
for (nm in c("categorical", "unadjusted")) {
  dfun <- if (nm == "categorical") design_categorical else design_unadjusted
  tabs <- lapply(setNames(cohort_order, cohort_order), function(k)
    fit_terms(prep[[k]]$dge, dfun(prep[[k]]$meta), "ballooning")$ballooning)
  a <- assemble(tabs)
  dt <- write_meta(a, file.path(out, sprintf("secondary_%s_meta_ballooning.tsv.gz", nm)))
  secondary[[nm]] <- list(
    n_meta_q05 = sum(dt$q_meta < 0.05), n_replicated = sum(dt$replicated),
    n_meta_q05_also_primary_q05 = length(intersect(dt[q_meta < 0.05, gene_id],
                                                   meta_obs$ballooning[q_meta < 0.05, gene_id])),
    n_replicated_also_primary_replicated = length(intersect(dt[replicated == TRUE, gene_id], rep_b$gene_id)),
    design_columns = lapply(setNames(cohort_order, cohort_order), function(k) colnames(dfun(prep[[k]]$meta))))
}

# ------------------------------------------------------- permutation null
set.seed(SEED)
strata_info <- list()
perm_index <- vector("list", n_perm)
strata_of <- lapply(setNames(cohort_order, cohort_order), function(k)
  paste(prep[[k]]$meta$fibrosis, prep[[k]]$meta$steatosis, sep = "_"))
for (k in cohort_order) {
  s <- strata_of[[k]]; tab <- table(s)
  strata_info[[k]] <- list(n_strata = length(tab), n_singleton_persons = sum(tab == 1L),
                           fraction_singleton = sum(tab == 1L) / length(s),
                           n_persons_in_strata_with_varying_ballooning =
                             sum(ave(prep[[k]]$meta$ballooning, s, FUN = function(z) length(unique(z))) > 1))
}
for (b in seq_len(n_perm)) {
  perm_index[[b]] <- lapply(setNames(cohort_order, cohort_order), function(k) {
    s <- strata_of[[k]]; idx <- seq_along(s)
    for (lv in unique(s)) {
      w <- which(s == lv)
      if (length(w) > 1L) idx[w] <- w[sample.int(length(w))]
    }
    idx
  })
}
run_perm <- function(b) {
  tabs <- lapply(setNames(cohort_order, cohort_order), function(k) {
    m <- copy(prep[[k]]$meta)
    m[, ballooning := ballooning[perm_index[[b]][[k]]]]
    fit_terms(prep[[k]]$dge, design_primary(m), "ballooning")$ballooning
  })
  a <- assemble(tabs)
  data.table(perm = b, n_replicated = sum(a$M$replicated), n_meta_q05 = sum(a$M$q_meta < 0.05))
}
t_perm <- Sys.time()
null <- rbindlist(mclapply(seq_len(n_perm), run_perm, mc.cores = n_cores, mc.preschedule = TRUE))
stopifnot(nrow(null) == n_perm)
fwrite(null, file.path(out, "permutation_null_counts.tsv"), sep = "\t")
perm_p <- (1 + sum(null$n_replicated >= N_obs)) / (1 + n_perm)
perm_p_q05 <- (1 + sum(null$n_meta_q05 >= sum(meta_obs$ballooning$q_meta < 0.05))) / (1 + n_perm)

answer <- if (!fib_pass) "not_interpretable_fibrosis_control_failed" else
  if (N_obs >= 1L && perm_p <= 0.05) "yes" else "no"

res_summary <- list(
  answer = answer,
  n_meta_family_genes = length(meta_genes),
  N_obs_replicated_ballooning = N_obs,
  n_ballooning_meta_q05 = sum(meta_obs$ballooning$q_meta < 0.05),
  replicated_counts = lapply(meta_obs, function(d) sum(d$replicated)),
  meta_q05_counts = lapply(meta_obs, function(d) sum(d$q_meta < 0.05)),
  within_cohort_q05 = lapply(setNames(cohort_order, cohort_order), function(k)
    lapply(obs[[k]], function(d) sum(d$q_within < 0.05))),
  genes_tested_per_cohort = lapply(prep, function(x) nrow(x$dge)),
  norm_method = lapply(prep, `[[`, "norm_method"),
  positive_control_fibrosis_pass = fib_pass,
  positive_control_steatosis_pass = ste_pass,
  permutation = list(n_perm = n_perm, seed = SEED, null_median = median(null$n_replicated),
                     null_q95 = unname(quantile(null$n_replicated, 0.95)), null_max = max(null$n_replicated),
                     perm_p = perm_p,
                     null_meta_q05_median = median(null$n_meta_q05),
                     null_meta_q05_q95 = unname(quantile(null$n_meta_q05, 0.95)),
                     perm_p_meta_q05 = perm_p_q05,
                     strata = strata_info,
                     minutes = as.numeric(difftime(Sys.time(), t_perm, units = "mins"))),
  secondary = secondary,
  minutes_total = as.numeric(difftime(Sys.time(), t_start, units = "mins")),
  candidate_not_adopted = TRUE)
writeLines(toJSON(res_summary, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(out, "summary.json"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
file.create(file.path(out, "COMPLETE"))
cat("done:", answer, "N_obs =", N_obs, "perm_p =", perm_p, "\n")
