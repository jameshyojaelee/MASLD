# DIG1 job 1 (exploratory; EXPLORATORY_PLAN.md questions A1-A3, B, C, D, E).
# Donor pseudobulks from V5 (GSE298719) and V3 (GSE256398). No cell-level data.
suppressPackageStartupMessages({library(data.table); library(edgeR)})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
OUT <- file.path(RUN, "out_pseudobulk_v2")
if (dir.exists(OUT)) stop("refusing to overwrite")
dir.create(OUT)
source(file.path(ROOT, "scripts/manuscript/program_observability_map/analysis_lib.R"))
V5 <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z")
J3 <- file.path(V5, "outputs/job3")
V3 <- file.path(ROOT, "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk")
HOT <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
BULK_PROG <- file.path(ROOT, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/per_program_all_arms.tsv")
H59 <- file.path(ROOT, "Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/results_primary_family.tsv")
BULK_DE <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
COLOC <- file.path(ROOT, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
DUCT <- "hotspot_hepatocytes_48f39dd4d817a10e"
ECM <- "hotspot_hepatocytes_f05c535ae5bbc0b9"
SEED <- 20261024
wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")

membership <- fread(file.path(HOT, "program_membership_v2.tsv"))
registry <- fread(file.path(HOT, "program_registry_v2.tsv"))
annotation <- fread(ANNOT, select = c("gene_id", "gene_name"))
annot_base <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
sym_of <- function(ids) annot_base$gene_name[match(sub("\\..*$", "", ids), annot_base$ensembl_base)]

# ---- cohorts and lineages ------------------------------------------------------------
dq <- fread(file.path(V5, "outputs/donor_qc.tsv"), colClasses = list(character = "id"))[included == TRUE][order(gsm)]
cellsF <- fread(file.path(J3, "cells_per_donor_set.tsv"))
detF <- fread(file.path(J3, "detection_by_set.tsv.gz"))
m3 <- fread(file.path(V3, "donor_metadata.tsv"))
n3 <- fread(file.path(V3, "nuclei_per_donor_lineage.tsv"))
donorsF <- data.table(donor = dq$gsm, disease = dq$disease, sex_f = dq$sex_f)
donors3 <- data.table(donor = m3$donor, disease = m3$disease, sex_f = as.integer(m3$sex == "female"))
COH <- list(GSE298719 = donorsF, GSE256398 = donors3)
lin <- data.table(
  cohort = c(rep("GSE298719", 5), rep("GSE256398", 3)),
  cell_type = c("hepatocytes", "cholangiocytes", "fibroblasts", "macrophages", "tcells",
                "hepatocytes", "cholangiocytes", "fibroblasts"),
  set = c("hep", "cholangiocytes", "fibroblasts", "macrophages", "tcells",
          "hepatocyte", "cholangiocyte", "stellate_fibroblast"))
lin[, file := ifelse(cohort == "GSE298719", file.path(J3, sprintf("pseudobulk_%s.tsv.gz", set)),
                     file.path(V3, sprintf("pseudobulk_%s.tsv.gz", set)))]
lin[, key := paste(cohort, cell_type, sep = ":")]

load_lineage <- function(i) {
  L <- lin[i]
  pb <- fread(L$file)
  meta <- copy(COH[[L$cohort]])
  if (L$cohort == "GSE298719") {
    meta[, n_cells := cellsF[set == L$set][match(meta$donor, gsm), n_cells]]
    keep <- detF$gene_id[detF[[L$set]] >= 0.01]
  } else {
    meta[, n_cells := n3[[L$set]][match(meta$donor, n3$donor)]]
    keep <- NULL
  }
  meta <- meta[n_cells >= 30]
  counts <- as.matrix(pb[, meta$donor, with = FALSE]); rownames(counts) <- pb$gene_id
  lc <- cpm(calcNormFactors(DGEList(counts)), log = TRUE, prior.count = 1)
  list(L = L, meta = meta, counts = counts, lc = lc, keep = keep)
}
LD <- setNames(lapply(seq_len(nrow(lin)), load_lineage), lin$key)
wr(rbindlist(lapply(LD, function(x) data.table(lkey = x$L$key, n_donors = nrow(x$meta),
  n_disease = sum(x$meta$disease), n_healthy = sum(x$meta$disease == 0), min_cells = min(x$meta$n_cells)))), "lineage_donors.tsv")

design_of <- function(meta, dis = meta$disease) cbind(1, disease = dis, sex_f = meta$sex_f, log10_n = log10(meta$n_cells))
hc3_mat <- function(Y, X, j = 2) {
  XtXi <- solve(crossprod(X)); C <- XtXi %*% t(X); B <- C %*% Y; E <- Y - X %*% B
  h <- rowSums((X %*% XtXi) * X); W <- E^2 / (1 - h)^2
  se <- sqrt(as.vector(C[j, , drop = FALSE]^2 %*% W)); df <- nrow(X) - ncol(X)
  data.table(beta = B[j, ], se = se, t = B[j, ] / se, df = df,
             ci_low = B[j, ] - qt(0.975, df) * se, ci_high = B[j, ] + qt(0.975, df) * se)
}
beta_only <- function(Y, X, j = 2) {
  if (qr(X)$rank < ncol(X)) return(rep(NA_real_, ncol(Y)))
  solve(crossprod(X), crossprod(X, Y))[j, ]
}
# Exact enumeration of cohort-wide disease relabelings (observed labeling included).
relabelings <- function(d) {
  cmb <- combn(nrow(d), sum(d$disease))
  apply(cmb, 2, function(ix) { v <- integer(nrow(d)); v[ix] <- 1L; v })
}
RL <- lapply(COH, relabelings)
wr(data.table(cohort = names(RL), n_relabelings = sapply(RL, ncol)), "relabeling_counts.tsv")
sub_lab <- function(x, v) v[match(x$meta$donor, COH[[x$L$cohort]]$donor)]

# ---- B: program scores and betas --------------------------------------------------------
score_lineage <- function(x) {
  reg <- registry[cell_type == x$L$cell_type]
  lc <- x$lc
  if (!is.null(x$keep)) lc <- lc[rownames(lc) %in% x$keep, , drop = FALSE]
  S <- collapse_symbols(lc, rownames(lc), annotation)
  sc <- score_programs(S, data.table(sample_id = x$meta$donor, dataset = x$L$cohort),
                       membership[program_uid %in% reg$program_uid], reg, 0.80)
  P <- sc$primary[, x$meta$donor, drop = FALSE]
  P[apply(P, 1, function(r) all(is.finite(r))), , drop = FALSE]
}
SC <- lapply(LD, score_lineage)
prog <- rbindlist(lapply(names(SC), function(k) {
  x <- LD[[k]]; Y <- t(SC[[k]])
  cbind(data.table(lkey = k, cohort = x$L$cohort, cell_type = x$L$cell_type, program_uid = colnames(Y),
                   n = nrow(Y)), hc3_mat(Y, design_of(x$meta)))
}))
prog <- merge(prog, registry[, .(program_uid, module, module_name)], by = "program_uid")
# Reproduction check against E6
e6 <- fread(file.path(J3, "E6_117_programs_descriptive.tsv"))[testable == TRUE]
chk <- merge(prog[cohort == "GSE298719", .(program_uid, beta)], e6[, .(program_uid, e6_beta = beta)], by = "program_uid")
wr(data.table(n_e6 = nrow(e6), n_matched = nrow(chk), max_abs_diff = max(abs(chk$beta - chk$e6_beta))), "check_E6_reproduction.tsv")
if (nrow(chk) != nrow(e6) || max(abs(chk$beta - chk$e6_beta)) > 1e-8) stop("E6 not reproduced")

# Discovery comparators
bp <- fread(BULK_PROG, select = c("feature_id", "new774_plus_GSE174478::beta_meta_fibrosis",
                                   "new774_plus_GSE174478::composition_pc_beta_meta_fibrosis",
                                   "new774_plus_GSE174478::q_value_fibrosis"))
setnames(bp, c("program_uid", "bulk_fib_beta", "bulk_fib_comp_beta", "bulk_fib_q"))
h59 <- fread(H59, select = c("cell_type", "module", "gls_beta", "gls_t"))
h59 <- merge(h59, registry[, .(cell_type, module, program_uid)], by = c("cell_type", "module"))[, .(program_uid, h59_beta = gls_beta, h59_t = gls_t)]
comp <- merge(merge(registry[, .(program_uid, cell_type, module_name)], bp, by = "program_uid", all.x = TRUE),
              h59, by = "program_uid", all.x = TRUE)
wide <- dcast(prog, program_uid ~ cohort, value.var = c("beta", "t"))
comp <- merge(comp, wide, by = "program_uid", all.x = TRUE)
wr(comp, "B_program_betas_all_sources.tsv")

spear <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 5) return(NA_real_); cor(a[ok], b[ok], method = "spearman") }
sign_agree <- function(a, b) { ok <- is.finite(a) & is.finite(b); mean(sign(a[ok]) == sign(b[ok])) }
# Null for one cohort's betas: relabel that cohort, keep comparators fixed.
null_betas <- function(coh) {
  ks <- names(SC)[lin$cohort[match(names(SC), lin$key)] == coh]
  R <- RL[[coh]]
  lapply(ks, function(k) {
    x <- LD[[k]]; Y <- t(SC[[k]])
    M <- vapply(seq_len(ncol(R)), function(r) beta_only(Y, design_of(x$meta, sub_lab(x, R[, r]))), numeric(ncol(Y)))
    if (is.null(dim(M))) M <- matrix(M, nrow = ncol(Y))
    rownames(M) <- colnames(Y); M
  }) |> setNames(ks)
}
NB <- lapply(names(COH), null_betas) |> setNames(names(COH))
conc_rows <- list()
for (coh in names(COH)) {
  other <- setdiff(names(COH), coh)
  comps <- list(bulk_fibrosis = comp$bulk_fib_beta, bulk_fibrosis_composition_adj = comp$bulk_fib_comp_beta,
                h59_singlecell_stage = comp$h59_beta)
  comps[[paste0("other_cohort_", other)]] <- comp[[paste0("beta_", other)]]
  B <- do.call(rbind, NB[[coh]])                       # programs x relabelings
  ct <- comp$cell_type[match(rownames(B), comp$program_uid)]
  obs_col <- which(apply(RL[[coh]], 2, function(v) all(v == COH[[coh]]$disease)))
  for (cn in names(comps)) {
    cv <- comps[[cn]][match(rownames(B), comp$program_uid)]
    for (grp in c("all", unique(ct))) {
      rows <- if (grp == "all") rep(TRUE, length(ct)) else ct == grp
      if (sum(rows & is.finite(cv)) < 5) next
      nullr <- apply(B[rows, , drop = FALSE], 2, spear, b = cv[rows])
      obs <- nullr[obs_col]
      conc_rows[[length(conc_rows) + 1]] <- data.table(
        cohort = coh, comparator = cn, lineage = grp, n_programs = sum(rows & is.finite(cv)),
        rho = obs, sign_agreement = sign_agree(B[rows, obs_col], cv[rows]),
        null_median = median(nullr, na.rm = TRUE), null_q95 = quantile(nullr, 0.95, na.rm = TRUE),
        p_upper_exact = mean(nullr >= obs - 1e-12, na.rm = TRUE), n_relabelings = sum(is.finite(nullr)))
    }
  }
}
wr(rbindlist(conc_rows), "B_concordance_exact_null.tsv")

# ---- gene-level DE (edgeR QL) ------------------------------------------------------------
de_lineage <- function(x) {
  y <- DGEList(x$counts)
  X <- model.matrix(~ disease + sex_f, data = x$meta)
  keep <- filterByExpr(y, design = X)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, X)
  fit <- glmQLFit(y, X, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = "disease"), n = Inf, sort.by = "none")$table
  data.table(lkey = x$L$key, gene_id = rownames(tt), symbol = sym_of(rownames(tt)), logFC = tt$logFC,
             logCPM = tt$logCPM, F = tt$F, P = tt$PValue, FDR = tt$FDR, signed_sqrtF = sign(tt$logFC) * sqrt(tt$F))
}
DE <- rbindlist(lapply(LD, de_lineage))
fwrite(DE, file.path(OUT, "gene_DE_edgeR_all_lineages.tsv.gz"), sep = "\t")
wr(DE[, .(n_genes = .N, fdr05 = sum(FDR < 0.05), fdr10 = sum(FDR < 0.10), up10 = sum(FDR < 0.10 & logFC > 0),
          down10 = sum(FDR < 0.10 & logFC < 0)), by = lkey], "gene_DE_summary.tsv")

# ---- A1: Ductular and ECM members, gene level -------------------------------------------
mem_d <- membership[program_uid == DUCT, .(symbol = mapped_symbol, weight = original_l1_weight)]
mem_e <- membership[program_uid == ECM, .(symbol = mapped_symbol, weight = original_l1_weight)]
a1 <- merge(DE[!is.na(symbol) & symbol %in% c(mem_d$symbol, mem_e$symbol)],
            rbind(mem_d[, program := "Ductular injury"], mem_e[, program := "Stromal ECM"]), by = "symbol", allow.cartesian = TRUE)
wr(a1[order(program, symbol, lkey)], "A1_member_gene_DE.tsv")

# ---- A2: cholangiocyte / stromal identity contamination in hepatocyte pseudobulk ----------
ident_genes <- function(coh, ref_ct, exclude, n_top = 30) {
  h <- LD[[paste0(coh, ":hepatocytes")]]; r <- LD[[paste0(coh, ":", ref_ct)]]
  dn <- intersect(h$meta[disease == 0, donor], r$meta[disease == 0, donor])
  basis <- "healthy"
  if (length(dn) < 3) { dn <- intersect(h$meta$donor, r$meta$donor); basis <- "all donors (fewer than 3 healthy pairs)" }
  g <- intersect(rownames(h$lc), rownames(r$lc))
  dif <- rowMeans(r$lc[g, dn, drop = FALSE] - h$lc[g, dn, drop = FALSE])
  ref_level <- rowMeans(r$lc[g, dn, drop = FALSE])
  tab <- data.table(gene_id = g, symbol = sym_of(g), diff = dif, ref_logcpm = ref_level)
  tab <- tab[!is.na(symbol) & !(symbol %in% exclude) & ref_logcpm >= 5][order(-diff)] |> head(n_top)
  tab[, `:=`(cohort = coh, reference = ref_ct, n_pairs = length(dn), basis = basis)]
  tab
}
a2_genes <- rbindlist(list(
  ident_genes("GSE298719", "cholangiocytes", mem_d$symbol),
  ident_genes("GSE256398", "cholangiocytes", mem_d$symbol),
  ident_genes("GSE298719", "fibroblasts", mem_e$symbol),
  ident_genes("GSE256398", "fibroblasts", mem_e$symbol)))
wr(a2_genes, "A2_identity_genes.tsv")
a2 <- list()
for (coh in names(COH)) {
  h <- LD[[paste0(coh, ":hepatocytes")]]
  for (ref in c("cholangiocytes", "fibroblasts")) {
    g <- a2_genes[cohort == coh & reference == ref, gene_id]
    contam <- colMeans(zscore_rows(h$lc[g, h$meta$donor, drop = FALSE]), na.rm = TRUE)
    prog_id <- if (ref == "cholangiocytes") DUCT else ECM
    sc <- if (prog_id %in% rownames(SC[[paste0(coh, ":hepatocytes")]])) SC[[paste0(coh, ":hepatocytes")]][prog_id, h$meta$donor] else rep(NA_real_, nrow(h$meta))
    X <- design_of(h$meta)
    c_fit <- hc3_mat(matrix(contam), X)
    nullc <- apply(RL[[coh]], 2, function(v) beta_only(matrix(contam), design_of(h$meta, sub_lab(h, v))))
    row <- data.table(cohort = coh, reference = ref, program = prog_id, contam_beta = c_fit$beta, contam_ci_low = c_fit$ci_low,
                      contam_ci_high = c_fit$ci_high, contam_p_exact = mean(abs(nullc) >= abs(c_fit$beta) - 1e-12, na.rm = TRUE),
                      rho_score_contam = if (all(is.finite(sc))) cor(sc, contam, method = "spearman") else NA_real_)
    if (all(is.finite(sc))) {
      f0 <- hc3_mat(matrix(sc), X); f1 <- hc3_mat(matrix(sc), cbind(X, contam = contam))
      row[, `:=`(score_beta = f0$beta, score_ci_low = f0$ci_low, score_ci_high = f0$ci_high,
                 score_beta_adj_contam = f1$beta, score_adj_ci_low = f1$ci_low, score_adj_ci_high = f1$ci_high)]
    }
    a2[[length(a2) + 1]] <- row
  }
}
wr(rbindlist(a2, fill = TRUE), "A2_contamination_control.tsv")

# ---- A3: lineage composition -------------------------------------------------------------
lcF <- fread(file.path(V5, "outputs/lineage_counts_per_donor.tsv"))
compF <- melt(lcF, id.vars = "gsm", variable.name = "lineage", value.name = "n")[, cohort := "GSE298719"]
setnames(compF, "gsm", "donor")
comp3 <- melt(n3, id.vars = "donor", variable.name = "lineage", value.name = "n")[, cohort := "GSE256398"]
cmp <- rbind(compF, comp3)
cmp[, frac := n / sum(n), by = .(cohort, donor)]
a3 <- cmp[, {
  d <- merge(data.table(donor = donor, frac = frac), COH[[.BY$cohort]], by = "donor")
  y <- qlogis(pmin(pmax(d$frac, 0.5 / 1e4), 1 - 0.5 / 1e4))
  X <- cbind(1, disease = d$disease, sex_f = d$sex_f)
  f <- hc3_mat(matrix(y), X)
  nl <- apply(RL[[.BY$cohort]], 2, function(v) beta_only(matrix(y), cbind(1, v[match(d$donor, COH[[.BY$cohort]]$donor)], d$sex_f)))
  list(median_frac_healthy = median(d$frac[d$disease == 0]), median_frac_disease = median(d$frac[d$disease == 1]),
       logit_beta = f$beta, ci_low = f$ci_low, ci_high = f$ci_high, p_exact = mean(abs(nl) >= abs(f$beta) - 1e-12, na.rm = TRUE))
}, by = .(cohort, lineage)]
wr(a3, "A3_composition.tsv")

# ---- C: GNMT, MAT1A, CYP2C19 and named genes ------------------------------------------------
named <- c("GNMT", "MAT1A", "CYP2C19", "FGF13", "MIR4435-2HG", "EPHA3", "ITGAV", "BICC1", "IGFBP7")
c_rows <- list()
for (k in names(LD)) {
  x <- LD[[k]]
  g <- rownames(x$lc)[sym_of(rownames(x$lc)) %in% named]
  if (!length(g)) next
  Y <- t(x$lc[g, x$meta$donor, drop = FALSE])
  f <- hc3_mat(Y, design_of(x$meta))
  N <- vapply(seq_len(ncol(RL[[x$L$cohort]])), function(r) beta_only(Y, design_of(x$meta, sub_lab(x, RL[[x$L$cohort]][, r]))), numeric(length(g)))
  if (is.null(dim(N))) N <- matrix(N, nrow = length(g))
  c_rows[[k]] <- cbind(data.table(lkey = k, gene_id = g, symbol = sym_of(g),
                                  mean_logcpm_healthy = rowMeans(x$lc[g, x$meta[disease == 0, donor], drop = FALSE])),
                       f, p_exact = rowMeans(abs(N) >= abs(f$beta) - 1e-12, na.rm = TRUE))
}
wr(rbindlist(c_rows), "C_named_genes_by_lineage.tsv")

# ---- D: lineage disease logFC versus the bulk disease transcriptome ----------------------------
bulk <- fread(BULK_DE)[, gene_id := sub("\\..*$", "", gene)]
bulk[, canonical := padj < 0.05 & abs(logFC) > 0.5]
d_rows <- list(); d_gene <- list()
for (k in names(LD)) {
  x <- LD[[k]]
  keep <- filterByExpr(DGEList(x$counts), design = model.matrix(~ disease + sex_f, data = x$meta))
  g <- intersect(rownames(x$lc)[keep], bulk$gene_id)
  Y <- t(x$lc[g, x$meta$donor, drop = FALSE])
  bl <- bulk$logFC[match(g, bulk$gene_id)]; canon <- bulk$canonical[match(g, bulk$gene_id)]
  R <- RL[[x$L$cohort]]
  stat <- function(b) c(rho = cor(b, bl, method = "spearman"), canon_sign = mean(sign(b[canon]) == sign(bl[canon])),
                        canon_up_pos = mean(b[canon & bl > 0] > 0), canon_down_neg = mean(b[canon & bl < 0] < 0))
  NM <- vapply(seq_len(ncol(R)), function(r) stat(beta_only(Y, design_of(x$meta, sub_lab(x, R[, r])))), numeric(4))
  obs_col <- which(apply(R, 2, function(v) all(v == COH[[x$L$cohort]]$disease)))
  o <- NM[, obs_col]
  d_rows[[k]] <- data.table(lkey = k, n_genes = length(g), n_canonical = sum(canon), stat = names(o), observed = o,
                            null_median = apply(NM, 1, median, na.rm = TRUE),
                            p_upper_exact = rowMeans(NM >= o - 1e-12, na.rm = TRUE))
  d_gene[[k]] <- data.table(lkey = k, gene_id = g, ols_beta = beta_only(Y, design_of(x$meta)))
}
wr(rbindlist(d_rows), "D_bulk_concordance_exact_null.tsv")
fwrite(rbindlist(d_gene), file.path(OUT, "D_gene_ols_beta.tsv.gz"), sep = "\t")

# ---- E: COLOC SuSiE-positive vs tested-negative genes -----------------------------------------
co <- fread(COLOC, select = c("ensembl", "susie_state_t12"))
co <- co[susie_state_t12 %in% c("positive", "tested_le_0.5")]
set.seed(SEED)
e_rows <- list()
for (k in unique(DE$lkey)) {
  tg <- DE[lkey == k & gene_id %in% co$ensembl]
  tg[, pos := co$susie_state_t12[match(gene_id, co$ensembl)] == "positive"]
  tg[, dec := cut(logCPM, quantile(logCPM, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
  tg[, a := abs(signed_sqrtF)]
  obs <- mean(tg$a[tg$pos]) - mean(tg$a[!tg$pos])
  nl <- replicate(10000, { p <- tg[, .(pos = sample(pos), a), by = dec]; mean(p$a[p$pos]) - mean(p$a[!p$pos]) })
  e_rows[[k]] <- data.table(lkey = k, n_positive = sum(tg$pos), n_negative = sum(!tg$pos),
                            mean_abs_t_positive = mean(tg$a[tg$pos]), mean_abs_t_negative = mean(tg$a[!tg$pos]),
                            diff = obs, p_two_sided_decile_perm = (1 + sum(abs(nl) >= abs(obs))) / 10001,
                            frac_fdr10_positive = mean(tg$FDR[tg$pos] < 0.10), frac_fdr10_negative = mean(tg$FDR[!tg$pos] < 0.10))
}
e <- rbindlist(e_rows)[, bh_q := p.adjust(p_two_sided_decile_perm, "BH")]
wr(e, "E_coloc_disease_change.tsv")
# Lineage of highest healthy expression (GSE298719, five lineages, genes filtered in each)
hl <- rbindlist(lapply(names(LD)[lin$cohort == "GSE298719"], function(k) {
  x <- LD[[k]]; data.table(gene_id = rownames(x$lc), cell_type = x$L$cell_type,
                           healthy_logcpm = rowMeans(x$lc[, x$meta[disease == 0, donor], drop = FALSE]))
}))
top <- hl[gene_id %in% co$ensembl][order(-healthy_logcpm)][, .SD[1], by = gene_id]
top[, state := co$susie_state_t12[match(gene_id, co$ensembl)]]
wr(dcast(top[, .N, by = .(state, cell_type)], cell_type ~ state, value.var = "N", fill = 0), "E_coloc_top_lineage_GSE298719.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
